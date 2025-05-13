# Load libraries
library(mgcv)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)  # for gradient_n_pal()

# Set seed for reproducibility
set.seed(42)

# Simulation parameters
n <- 500
t <- seq(0, 99, length.out = n)
T0 <- 0
k_T <- -0.5
k_D <- 2.0
sigma <- 0.5

# Temperature variation: 10-period sinusoidal
T_t <- T0 + 5 * sin(2 * pi * 5 * t / max(t)) + 0.02 *t 

# Quadratic, centered temperature effect
temp_effect <- k_T * (T_t - T0)^2

# Drift and noise
drift_effect <- k_D * sqrt(t)
noise <- rnorm(n, 0, sigma)

# Observed signal
R_t <- temp_effect + drift_effect + noise

# Data frame
df <- data.frame(
  R = R_t,
  temp = T_t,
  time = t,
  temp_effect = temp_effect,
  drift_effect = drift_effect,
  noise = noise
)

# Fit GAM model
model <- gam(R ~ s(temp, bs = "cr", k=10) + s(time, pc=0, bs = "cr", k = 10 ), data = df)
# model <- gam(R ~ s(temp, bs = "ps", k=10) + s(time, pc=0, bs = "ps", k = 10 ), data = df)

# Create prediction sequences
temp_seq <- seq(min(df$temp), max(df$temp), length.out = 200)
time_seq <- seq(min(df$time), max(df$time), length.out = 200)

# --- Predict temperature effect ---
df_temp_pred <- data.frame(temp = temp_seq, time = mean(df$time))
temp_pred_obj <- predict(model, newdata = df_temp_pred, type = "terms", terms = "s(temp)", se.fit = TRUE)
temp_pred <- temp_pred_obj$fit[, "s(temp)"]
temp_se   <- temp_pred_obj$se.fit[, "s(temp)"]

# --- Predict time effect ---
df_time_pred <- data.frame(temp = mean(df$temp), time = time_seq)
time_pred_obj <- predict(model, newdata = df_time_pred, type = "terms", terms = "s(time)", se.fit = TRUE)
time_pred <- time_pred_obj$fit[, "s(time)"]
time_se <- time_pred_obj$se.fit[, "s(time)"]

# Define IEEE-style theme
ieee_theme <- theme_minimal(base_size = 10, base_family = "sans") +
  theme(
    plot.title = element_text(size = 11, face = "bold"),
    legend.position = "right",
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_blank()
  )

# --- Plot 1: Raw component decomposition ---
df_long <- df %>%
  select(time, R, temp_effect, drift_effect, noise) %>%
  pivot_longer(cols = -time, names_to = "component", values_to = "value")

p1 <- ggplot(df_long, aes(x = time, y = value, color = component)) +
  geom_line(size = 1) +
  scale_color_manual(
    values = c("R" = "black", "temp_effect" = "forestgreen", "drift_effect" = "blue", "noise" = "red"),
    labels = c("R" = "Measurement error", "temp_effect" = "Temperature Effect",
               "drift_effect" = "Instrumental drift", "noise" = "Noise")
  ) +
  labs(x = "Time (days)", y = "Measurement error (g)", color = "Component") +
  ieee_theme

# --- Plot 2: Temp effect (with error band) ---
true_temp_seq <- k_T * (temp_seq - T0)^2
true_temp_seq <- true_temp_seq - mean(true_temp_seq)

df_temp_effect <- data.frame(
  temp = temp_seq,
  True = true_temp_seq,
  Estimated = temp_pred,
  Upper = temp_pred + 1.96 * temp_se,
  Lower = temp_pred - 1.96 * temp_se
)

p2 <- ggplot(df_temp_effect, aes(x = temp)) +
  geom_line(aes(y = True, color = "True"), size = 1, linetype = "dashed") +
  geom_line(aes(y = Estimated, color = "Estimated"), size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "Estimated"), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  scale_fill_manual(values = c("Estimated" = "blue")) +
  labs(x = "Temperature (\u00B0C)", y = "Effect (g)", color = "Effect Type", fill = "Confidence") +
  ieee_theme

# --- Plot 3: Time effect (with error band) ---
true_time_seq <- k_D * sqrt(time_seq)

df_time_effect <- data.frame(
  time = time_seq,
  True = true_time_seq,
  Estimated = time_pred,
  Upper = time_pred + 1.96 * time_se,
  Lower = time_pred - 1.96 * time_se
)

p3 <- ggplot(df_time_effect, aes(x = time)) +
  geom_line(aes(y = True, color = "True"), size = 1, linetype = "dashed") +
  geom_line(aes(y = Estimated, color = "Estimated"), size = 1) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper, fill = "Estimated"), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  scale_fill_manual(values = c("Estimated" = "blue")) +
  labs(x = "Time (days)", y = "Effect (g)", color = "Effect Type", fill = "Confidence") +
  ieee_theme

# --- Plot 4: Temp effect derivative ---
temp_deriv_true <- diff(true_temp_seq) / diff(temp_seq)
temp_deriv_est  <- diff(temp_pred) / diff(temp_seq)

df_temp_deriv <- data.frame(
  temp = temp_seq[-1],
  Derivative = c(temp_deriv_true, temp_deriv_est),
  Type = rep(c("True", "Estimated"), each = length(temp_deriv_true))
)

p4 <- ggplot(df_temp_deriv, aes(x = temp, y = Derivative, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  scale_linetype_manual(values = c("True" = "dashed", "Estimated" = "solid")) +
  labs(x = "Temperature (\u00B0C)", y = "First Derivative (g/\u00B0C)", color = "Effect Type", linetype = "Effect Type") +
  ieee_theme

# --- Plot 5: Time effect derivative ---
time_deriv_true <- diff(true_time_seq) / diff(time_seq)
time_deriv_est  <- diff(time_pred) / diff(time_seq)

df_time_deriv <- data.frame(
  time = time_seq[-1],
  Derivative = c(time_deriv_true, time_deriv_est),
  Type = rep(c("True", "Estimated"), each = length(time_deriv_true))
)

p5 <- ggplot(df_time_deriv, aes(x = time, y = Derivative, color = Type, linetype = Type)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("True" = "black", "Estimated" = "blue")) +
  scale_linetype_manual(values = c("True" = "dashed", "Estimated" = "solid")) +
  labs(x = "Time (days)", y = "First Derivative (g/day)", color = "Effect Type", linetype = "Effect Type") +
  ieee_theme

# --- Plot 6: Measurement error vs temperature, color-coded by time ---
p6 <- ggplot(df, aes(x = temp, y = R, color = time)) +
  geom_point(size = 0.7) +
  scale_color_gradientn(colors = c("#00B945", "#0C5DA5"), name = "Time") +
  labs(x = "Temperature (\u00B0C)", y = "Measurement Error (g)") +
  ieee_theme

# --- Plot 7: Temperature vs Time ---
p7 <- ggplot(df, aes(x = time, y = temp)) +
  geom_line(color = "blue", size = 1) +
  labs(x = "Time (days)", y = "Temperature (\u00B0C)") +
  ieee_theme

# --- Save individual plots ---
plot_list <- list(p1 = p1, p2 = p2, p3 = p3, p4 = p4, p5 = p5, p6 = p6, p7 = p7)
plot_names <- names(plot_list)

# Get the path of the currently running script (robust fallback)
get_script_path <- function() {
  # Works when sourced via Rscript or using RStudio
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else if (!is.null(sys.frames()[[1]]$ofile)) {
    return(normalizePath(sys.frames()[[1]]$ofile))
  } else {
    stop("Cannot determine script location. Please set output directory manually.")
  }
}

script_path <- get_script_path()
script_dir <- dirname(script_path)
output_dir <- file.path(script_dir, "figs R")

# Create directory if it doesn't exist
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Replace plot_names with your custom file names
file_names <- c(
  "syn-err-decomp", "syn-est-temp", "syn-est-time",
  "syn-est-temp-dt", "syn-est-time-dt", "syn-err-vs-temp", "syn-temp"
)

# Save each plot with the correct filename
for (i in seq_along(plot_list)) {
  ggsave(
    filename = file.path(output_dir, paste0(file_names[i], ".png")),
    plot = plot_list[[i]],
    width = 4.5, height = 3.5, dpi = 600, device = "png"
  )
}
