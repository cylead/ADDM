library(ggplot2)
library(scales)

# ---- Setup ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ---- Load saved data ----
load("./saved_data/R/gam_effects_data.RData")

# ---- IEEE Transactions Theme ----
ieee_theme <- theme_minimal() + 
  theme(
    legend.position = "bottom",
    legend.key.size = unit(0.4, "cm"),
    legend.text = element_text(size=7),
    axis.title = element_text(size=8),
    axis.text = element_text(size=7),
    axis.text.x = element_text(angle=45, hjust=1),
    plot.margin = margin(2, 2, 2, 2, "mm"),  # Reduced margins
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(size=0.2)
  )

# ---- Compute means and std ----
T_mean <- rowMeans(T_mat)
T_sd <- apply(T_mat, 1, sd)

t_mean <- rowMeans(t_mat)
t_sd <- apply(t_mat, 1, sd)

RH_mean <- rowMeans(RH_mat)
RH_sd <- apply(RH_mat, 1, sd)

# ---- Create output directory ----
if(!dir.exists("./figs/GAM R")) dir.create("./figs/GAM R", recursive = TRUE)

# ----------------------------
# ---- Combined Effect Plots --
# ----------------------------

# ---- Combined Time Effect ----
t0_global <- t0
t_test_timestamps_global <- t0_global + t_range

df_time <- data.frame(
  x = rep(t_test_timestamps_global, length(labels_used)),
  y = as.vector(t_mat)
)

p_comb_time <- ggplot() +
  geom_line(data = df_time, aes(x = x, y = y, group=rep(labels_used, each=100)), color="grey70", alpha=0.5) +
  geom_line(aes(x=t_test_timestamps_global, y=t_mean), color="black") +
  geom_ribbon(aes(x=t_test_timestamps_global, ymin=t_mean - 2*t_sd, ymax=t_mean + 2*t_sd), fill="black", alpha=0.2) +
  labs(x="Time", y="Instrumental drift (ppm)") +
  scale_x_datetime(date_breaks = "3 months", date_labels = "%Y-%m") +
  scale_y_continuous(breaks = pretty, n.breaks = 6) +
  ieee_theme

ggsave("./figs/GAM R/Combined_Time.png", p_comb_time, width=3.5, height=4, dpi=600, units = "in")

# ---- Combined Temperature Effect ----
df_temp <- data.frame(
  x = rep(T_range, length(labels_used)),
  y = as.vector(T_mat)
)

p_comb_temp <- ggplot() +
  geom_line(data = df_temp, aes(x = x, y = y, group=rep(labels_used, each=100)), color="grey70", alpha=0.5) +
  geom_line(aes(x=T_range, y=T_mean), color="black") +
  geom_ribbon(aes(x=T_range, ymin=T_mean - 2*T_sd, ymax=T_mean + 2*T_sd), fill="black", alpha=0.2) +
  labs(x="Temperature (\u00B0C)", y="Effect (ppm)") +
  ieee_theme

ggsave("./figs/GAM R/Combined_Temperature.png", p_comb_temp, width=3.5, height=4, dpi=600, units = "in")

# ---- Combined RH Effect ----
df_rh <- data.frame(
  x = rep(RH_range, length(labels_used)),
  y = as.vector(RH_mat)
)

p_comb_rh <- ggplot() +
  geom_line(data = df_rh, aes(x = x, y = y, group=rep(labels_used, each=100)), color="grey70", alpha=0.5) +
  geom_line(aes(x=RH_range, y=RH_mean), color="black") +
  geom_ribbon(aes(x=RH_range, ymin=RH_mean - 2*RH_sd, ymax=RH_mean + 2*RH_sd), fill="black", alpha=0.2) +
  labs(x="Relative Humidity (%)", y="Effect (ppm)") +
  ieee_theme

ggsave("./figs/GAM R/Combined_RH.png", p_comb_rh, width=3.5, height=4, dpi=600, units = "in")

cat("Saved combined plots WITHOUT curve labels into ./figs/GAM R/\n")
