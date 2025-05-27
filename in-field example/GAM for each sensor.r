library(mgcv)
library(ggplot2)
library(dplyr)
library(readr)
library(scales)
library(lubridate)

# ---- Setup ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Load data ----
data <- read_csv("./saved_data/R/gam_data_for_R.csv")

# ---- Create output directory if it doesn't exist ----
if(!dir.exists("./figs/GAM R")) dir.create("./figs/GAM R", recursive = TRUE)

# ---- IEEE Transactions Theme ----
ieee_theme <- theme_minimal() +
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size=7),
    axis.title = element_text(size=8),
    axis.text = element_text(size=7),
    axis.text.x = element_text(angle=45, hjust=1),
    plot.margin = unit(c(0.3, 0.3, 0.3, 0.3), "cm"),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size=0.2)
  )

t0 <- min(data$timestamp)   # <<<< Set once for all labels!

# ---- Define trouble list ----
trouble_list <- c("1076", "1141", "1173", "1179", "1278")

# ---- Initialize storage for combined plots ----
T_effects_all <- list()
t_effects_all <- list()
RH_effects_all <- list()
labels_used <- c()

# ---- Loop over all labels ----
for(label in setdiff(unique(data$label), trouble_list)) {
  
  cat("Processing label:", label, "\n")
  
  df <- filter(data, label == !!label)
  
  if(nrow(df) < 50) next  # Skip if too few data points
  
  labels_used <- c(labels_used, label)
  
  # ---- Time conversion ----
  t_timestamps <- t0 + df$t
  
  # ---- Fit GAM ----
  # gam_model <- gam(err ~ s(T) + s(t, k=30, pc=0) + s(RH), data=df)
  gam_model <- gam(err ~ s(T, bs = "cr", k=20) + s(t, bs='cr', k=200, pc=0) + 
                                        s(RH, bs='ad', k=30), data=df)
  
  # ---- Define ranges ----
  T_range <- seq(min(df$T), max(df$T), length.out=100)
  t_range <- seq(min(df$t), max(df$t), length.out=100)
  RH_range <- seq(min(df$RH), max(df$RH), length.out=100)
  t_test_timestamps <- t0 + t_range
  
  # ---- Partial Effects ----
  T_test <- data.frame(T=T_range, t=mean(df$t), RH=mean(df$RH))
  t_test <- data.frame(T=mean(df$T), t=t_range, RH=mean(df$RH))
  RH_test <- data.frame(T=mean(df$T), t=mean(df$t), RH=RH_range)
  
  T_terms <- predict(gam_model, newdata=T_test, type="terms", se.fit=TRUE)
  t_terms <- predict(gam_model, newdata=t_test, type="terms", se.fit=TRUE)
  RH_terms <- predict(gam_model, newdata=RH_test, type="terms", se.fit=TRUE)
  
  T_effects_all[[label]] <- T_terms$fit[,1]
  t_effects_all[[label]] <- t_terms$fit[,2]
  RH_effects_all[[label]] <- RH_terms$fit[,3]
  
  # ---- Save Individual Plots ----
  plots <- list(
    Time = list(x = t_test_timestamps, effect = t_terms$fit[,2], lower = t_terms$fit[,2] - 2*t_terms$se.fit[,2], upper = t_terms$fit[,2] + 2*t_terms$se.fit[,2], xlab="Time", ylab="Instrumental drift (ppm)"),
    Temperature = list(x = T_range, effect = T_terms$fit[,1], lower = T_terms$fit[,1] - 2*T_terms$se.fit[,1], upper = T_terms$fit[,1] + 2*T_terms$se.fit[,1], xlab="Temperature (\u00B0C)", ylab="Effect (ppm)"),
    RH = list(x = RH_range, effect = RH_terms$fit[,3], lower = RH_terms$fit[,3] - 2*RH_terms$se.fit[,3], upper = RH_terms$fit[,3] + 2*RH_terms$se.fit[,3], xlab="Relative Humidity (%)", ylab="Effect (ppm)")
  )
  
  for(pn in names(plots)){
    p <- plots[[pn]]
    plt <- ggplot() +
      geom_ribbon(aes(x=p$x, ymin=p$lower, ymax=p$upper, fill="95% CI"), alpha=0.2) +
      geom_line(aes(x=p$x, y=p$effect, color="Estimated effect")) +
      geom_hline(aes(yintercept=0, linetype="y=0 reference"), color="red") +
      labs(x=p$xlab, y=p$ylab) +
      scale_color_manual("", values=c("Estimated effect"="black")) +
      scale_fill_manual("", values=c("95% CI"="black")) +
      scale_linetype_manual("", values=c("y=0 reference"="dotted")) +
      ieee_theme
    
    ggsave(paste0("./figs/GAM R/", label, "_", pn, ".png"), plt, width=3.5, height=2.5, dpi=600)
  }
  
  # ---- Plot Residuals ----
  df$residuals <- residuals(gam_model)
  
  p_res <- ggplot(df, aes(x=t_timestamps, y=residuals)) +
    geom_point(aes(color="Residuals"), alpha=0.7, size=0.1) +
    geom_hline(aes(yintercept=0, linetype="y=0 reference"), color="red") +
    labs(x="Time", y="Residuals (ppm)") +
    ylim(-100, 100) +
    scale_color_manual("", values=c("Residuals"="black")) +
    scale_linetype_manual("", values=c("y=0 reference"="dotted")) +
    scale_x_datetime(date_breaks = "3 months", date_labels = "%Y-%m") +
    ieee_theme
  
  ggsave(paste0("./figs/GAM R/", label, "_Residuals.png"), p_res, width=3.5, height=2.5, dpi=600)
}

# ---------------------------
# ---- Combined Effect Plots ----
# ---------------------------

# --- Convert to matrices ---
T_mat <- do.call(cbind, T_effects_all)
t_mat <- do.call(cbind, t_effects_all)
RH_mat <- do.call(cbind, RH_effects_all)

# --- Compute means and std ---
T_mean <- rowMeans(T_mat)
T_sd <- apply(T_mat, 1, sd)

t_mean <- rowMeans(t_mat)
t_sd <- apply(t_mat, 1, sd)

RH_mean <- rowMeans(RH_mat)
RH_sd <- apply(RH_mat, 1, sd)

save(
  T_range, t_range, RH_range,
  T_mat, t_mat, RH_mat,
  labels_used,
  t0,
  file = "./saved_data/R/gam_effects_data.RData"
)

cat("Saved all effect data to ./saved_data/R/gam_effects_data.RData\n")
