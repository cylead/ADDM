rm(list = ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

library(readr)
library(dplyr)
library(rlang)

##----- load data -----
data <- read_csv("./saved_data/R/gam_augmented_all_labels.csv",
                 show_col_types = FALSE)
out_dir <- "./figs R/revision 1/improved calibration"
# Columns assumed present in 'data':
# lp8, effect_T, effect_RH, effect_T_avg, effect_RH_avg,
# ref, err, ts_round, label, temperature

label_list <- c("1187", "1215", "1029", "1045", "1312", "1288",
                "1022", "1097", "1098", "1094", "1103", "1104", "1055")

##----- helper: RMSE -----
rmse <- function(x) sqrt(mean(x^2, na.rm = TRUE))

##----- generic calibration function for one scheme -----
calibrate_scheme <- function(df, value_col,
                             use_ref_target = FALSE,
                             target_value = 400) {
  v <- sym(value_col)
  
  # find, per week, the row with the lowest value_col
  mins <- df %>%
    group_by(week) %>%
    slice_min(order_by = !!v, with_ties = FALSE) %>%
    ungroup() %>%
    transmute(
      week,
      baseline     = !!v,
      ref_baseline = ref
    )
  
  df_cal <- df %>%
    left_join(mins, by = "week") %>%
    mutate(
      offset = if (use_ref_target) {
        # not used here, always calibrate to fixed target
        ref_baseline - baseline
      } else {
        target_value - baseline
      },
      calibrated = !!v + offset,
      error_cal  = calibrated - ref
    )
  
  rmse(df_cal$error_cal)
}

##----- function to compute 7 RMSE schemes for a single label -----
# Schemes:
# 1) lp8_to_400: lp8 → weekly minimum becomes 400
# 2) lp8_minusT_to_400: (lp8 - effect_T) → min = 400
# 3) lp8_minusRH_to_400: (lp8 - effect_RH) → min = 400
# 4) lp8_minusT_minusRH_to_400: (lp8 - effect_T - effect_RH) → min = 400
# 5) lp8_minusT_avg_to_400: (lp8 - effect_T_avg) → min = 400
# 6) lp8_minusRH_avg_to_400: (lp8 - effect_RH_avg) → min = 400
# 7) lp8_minusT_avg_RH_avg_to_400: (lp8 - effect_T_avg - effect_RH_avg) → min = 400
compute_rmse_for_label <- function(df_all, lbl) {
  df_lbl <- df_all %>%
    filter(label == lbl) %>%
    # define week index: non-overlapping 7-day blocks
    mutate(
      week = floor(as.numeric(ts_round - min(ts_round, na.rm = TRUE)) /
                     (7 * 24 * 3600))
    ) %>%
    # corrected signals
    mutate(
      lp8_Tcorr           = lp8 - effect_T,
      lp8_RHcorr          = lp8 - effect_RH,
      lp8_TRHcorr         = lp8 - effect_T - effect_RH,
      lp8_Tavg_corr       = lp8 - effect_T_avg,
      lp8_RHavg_corr      = lp8 - effect_RH_avg,
      lp8_Tavg_RHavg_corr = lp8 - effect_T_avg - effect_RH_avg
    )
  
  ##----- seven schemes (all to fixed target 400) -----
  rmse_scheme1 <- calibrate_scheme(df_lbl, "lp8",
                                   use_ref_target = FALSE, target_value = 400)
  rmse_scheme2 <- calibrate_scheme(df_lbl, "lp8_Tcorr",
                                   use_ref_target = FALSE, target_value = 400)
  rmse_scheme3 <- calibrate_scheme(df_lbl, "lp8_RHcorr",
                                   use_ref_target = FALSE, target_value = 400)
  rmse_scheme4 <- calibrate_scheme(df_lbl, "lp8_TRHcorr",
                                   use_ref_target = FALSE, target_value = 400)
  rmse_scheme5 <- calibrate_scheme(df_lbl, "lp8_Tavg_corr",
                                   use_ref_target = FALSE, target_value = 400)
  rmse_scheme6 <- calibrate_scheme(df_lbl, "lp8_RHavg_corr",
                                   use_ref_target = FALSE, target_value = 400)
  rmse_scheme7 <- calibrate_scheme(df_lbl, "lp8_Tavg_RHavg_corr",
                                   use_ref_target = FALSE, target_value = 400)
  
  ##----- RMSE of original err column (optional) -----
  rmse_err <- rmse(df_lbl$err)
  
  ##----- collect -----
  out <- data.frame(
    label                           = lbl,
    lp8_to_400                      = rmse_scheme1,
    lp8_minusT_to_400               = rmse_scheme2,
    lp8_minusRH_to_400              = rmse_scheme3,
    lp8_minusT_minusRH_to_400       = rmse_scheme4,
    lp8_minusT_avg_to_400           = rmse_scheme5,
    lp8_minusRH_avg_to_400          = rmse_scheme6,
    lp8_minusT_avg_RH_avg_to_400    = rmse_scheme7,
    original_err                    = rmse_err
  )
  
  out
}

##---------------- ALL TEMPERATURES ----------------##
rmse_by_label_list <- lapply(label_list, function(lbl) {
  compute_rmse_for_label(data, lbl)
})

rmse_by_label <- do.call(rbind, rmse_by_label_list)

# compute column means (excluding label column)
mean_row <- rmse_by_label %>%
  select(-label) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  mutate(label = "MEAN") %>%
  select(label, everything())

rmse_by_label <- bind_rows(rmse_by_label, mean_row)

##---------------- PLOTTING (IEEE STYLE, PNG) ----------------##

# Use the 7 scheme columns (2:8) – excluding 'original_err'
# rmse_base_all     <- rmse_by_label[rmse_by_label$label != "MEAN", 2:8]

rmse_base_all <- rmse_by_label %>%
  filter(label != "MEAN") %>%
  select(
    lp8_to_400,
    lp8_minusT_minusRH_to_400,
    lp8_minusT_avg_RH_avg_to_400
  )

mean_row <- rmse_base_all %>%
  summarise(across(everything(), ~round(mean(.x, na.rm = TRUE), 2))) %>%
  mutate(stat = "mean", .before = 1)

median_row <- rmse_base_all %>%
  summarise(across(everything(), ~round(median(.x, na.rm = TRUE), 2))) %>%
  mutate(stat = "median", .before = 1)
rmse_stats_table <- bind_rows(mean_row, median_row)

# Save to CSV in the SAME directory as image output
write_csv(
  rmse_stats_table,
  file.path(out_dir, "rmse_stats_mean_median.csv")
)

# Short scheme names for x-axis (you can tweak as you like)
# scheme_names <- c(
#   "ABC",        # lp8_to_400
#   "T-ABC",      # lp8_minusT_to_400
#   "RH-ABC",     # lp8_minusRH_to_400
#   "E-ABC",    # lp8_minusT_minusRH_to_400
#   "Tavg-ABC",   # lp8_minusT_avg_to_400
#   "RHavg-ABC",  # lp8_minusRH_avg_to_400
#   "ET-ABC"  # lp8_minusT_avg_RH_avg_to_400
# )

scheme_names <- c(
  "ABC",        # lp8_to_400
  "E-ABC",    # lp8_minusT_minusRH_to_400
  "EE-ABC"  # lp8_minusT_avg_RH_avg_to_400
)
##----- create output directory -----
out_dir <- "./figs R/revision 1/improved calibration"

##========================================================##
##  PNG 1: ALL TEMPERATURES
##========================================================##

png(
  filename = file.path(out_dir, "rmse_box.png"),
  width  = 3,   # in inches (typical IEEE figure width per subfigure)
  height = 2,
  units  = "in",
  res    = 600    # high resolution for publication
)

par(
  family = "serif",
  mar = c(2.5, 3.0, 0.5, 0.5),
  oma = c(0,0,0,0),
  mgp = c(1.8, 0.5, 0)
)
boxplot(
  rmse_base_all,
  ylim  = c(10, 35),       # fixed y-range for comparability
  names = scheme_names,    # custom x labels
  las   = 1,               # rotate labels if needed
  main  = "",
  xlab  = "",
  ylab  = "RMSE",
  col   = "white",
  border = "black",
  lwd   = 1,
  cex.axis = 0.8,
  cex.lab  = 1
)

dev.off()

# Select only the three schemes and keep labels
rmse_all_three <- rmse_by_label %>%
  filter(label != "MEAN") %>%
  select(
    label,
    lp8_to_400,
    lp8_minusT_minusRH_to_400,
    lp8_minusT_avg_RH_avg_to_400
  ) %>%
  # round ALL numeric columns to 3 decimals
  mutate(across(where(is.numeric), ~ round(.x, 3)))

# Compute mean row
mean_row <- rmse_all_three %>%
  summarise(across(-label, ~ round(mean(.x, na.rm = TRUE), 3))) %>%
  mutate(label = "MEAN") %>%
  select(label, everything())

# Compute median row
median_row <- rmse_all_three %>%
  summarise(across(-label, ~ round(median(.x, na.rm = TRUE), 3))) %>%
  mutate(label = "MEDIAN") %>%
  select(label, everything())

# Combine all rows
rmse_stats_all <- bind_rows(
  rmse_all_three,
  mean_row,
  median_row
)
rmse_stats_all <- rmse_stats_all %>%
  mutate(across(
    where(is.numeric),
    ~ format(round(.x, 3), nsmall = 3, trim = TRUE)
  ))

write_csv(
  rmse_stats_all,
  file.path(out_dir, "rmse_stats_all_labels_with_mean_median.csv")
)

