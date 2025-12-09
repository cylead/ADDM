rm(list = ls())

library(mgcv)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)
library(purrr)
library(data.table)

# --------- ACF for irregular (duplicate/gappy) time index ----------
# Computes autocorrelation at lag h as correlation between *all observations*
# whose time difference is h days. Duplicates are handled naturally, gaps
# simply contribute no pairs.
acf_irregular <- function(residuals, t, lag.max = 10, unbiased = TRUE) {
  # residuals: numeric residual vector
  # t: integer (or coercible) time index, e.g., days (can have duplicates, NAs, gaps)
  # lag.max: maximum lag
  # unbiased = FALSE:
  #   acf(h) = S_h / S_0, where S_h = sum_{pairs at lag h} (x_i - mu)(x_j - mu),
  #   S_0 = sum_i (x_i - mu)^2
  # unbiased = TRUE:
  #   acf(h) = (S_h / N_h) / (S_0 / N), i.e., covariance with lag-specific pair count N_h
  #   divided by overall variance.
  
  stopifnot(length(residuals) == length(t))
  
  DT <- data.table(t = as.integer(t), e = as.numeric(residuals))
  DT <- DT[!is.na(t) & !is.na(e)]
  if (nrow(DT) == 0L) stop("No non-missing residuals/t values.")
  
  mu <- mean(DT$e)
  # Aggregate by day once: sums and counts
  agg <- DT[, .(sum_e = sum(e), n = .N), by = t]
  setkey(agg, t)
  
  # Basic totals for normalization
  N   <- nrow(DT)
  S0  <- sum((DT$e - mu)^2)  # denominator for acf if unbiased = FALSE
  var_unbiased <- S0 / N     # for unbiased-style normalization
  
  # Prepare result container
  out <- data.table(
    lag   = 0:lag.max,
    acf   = NA_real_,
    pairs = NA_real_  # number of contributing pairs N_h
  )
  
  # Lag 0: define ACF(0) = 1 by convention
  out[lag == 0L, `:=`(acf = 1, pairs = as.numeric(N))]
  
  # For lags > 0, compute efficiently via aggregated sums
  for (h in seq_len(lag.max)) {
    # Make a shifted copy so that we can inner-join on matching base day t
    b <- copy(agg)
    b[, t := t - h]  # rows in b come from original day (t+h); shift to align with 'agg' day t
    
    # Join: keeps only days where both t and t+h exist
    m <- agg[b, on = "t", nomatch = 0L]
    if (nrow(m) == 0L) {
      out[lag == h, `:=`(acf = NA_real_, pairs = 0)]
      next
    }
    
    # Centered group sums for day d and day d+h
    sx <- m$sum_e - m$n * mu
    sy <- m$i.sum_e - m$i.n * mu
    
    # Numerator and pair counts
    S_h <- sum(sx * sy)
    N_h <- sum(m$n * m$i.n)
    
    if (!unbiased) {
      acf_h <- if (S0 == 0) NA_real_ else S_h / S0
    } else {
      # (covariance at lag h) / (overall variance)
      acf_h <- if (N_h == 0 || var_unbiased == 0) NA_real_ else (S_h / N_h) / var_unbiased
    }
    
    out[lag == h, `:=`(acf = acf_h, pairs = as.numeric(N_h))]
  }
  
  out[]
}

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# ensure output directory exists
plot_dir <- "./figs R/revision 1/basis gam for all"

# Set the origin date corresponding to t_round (days since origin)
origin_date <- as.Date("2017-07-01")

data <- read_csv("./saved_data/R/gam_data_outlier_removed.csv", show_col_types = FALSE)

label_list <- c(
  "1187", "1215", "1029", "1045", "1312", "1288", "1022",
  "1097", "1098", "1094", "1103", "1104", "1055"
)
# label_list <- c("1097", "1215")

# containers
all_augmented <- list()
model_stats   <- list()
models        <- list()

# IEEE-style figure size (single-column) and resolution
ieee_width  <- 3.5  # inches (single column)
ieee_height <- 2.5  # inches
ieee_dpi    <- 600  # high resolution for print

# Common base theme for IEEE-style figures (smaller text, suitable for 3.5" width)
base_theme_ieee <- theme_bw(base_size = 8)

for (lbl in label_list) {
  message("Processing label: ", lbl)
  
  # Subset data: use existing t_round from the dataset, just sort by timestamp
  df_lbl <- dplyr::filter(data, label == lbl) %>%
    arrange(timestamp)
  
  # -------------------------------------------------------------------
  # Fit GAM (no heteroscedasticity structure)
  # -------------------------------------------------------------------
  gm <- bam(
    err ~ s(T, k = 10, bs = "cr") +
      s(t_round, k = 200, bs = "cr", pc = 0) +
      s(RH, k = 10, bs = "cr"),
    data   = df_lbl, method = "fREML", discrete = TRUE
  )
  
  # store model
  models[[lbl]] <- gm
  
  # -------------------------------------------------------------------
  # Add effects + residuals to df_lbl
  # -------------------------------------------------------------------
  term_mat <- predict(gm, newdata = df_lbl, type = "terms")
  term_df  <- as.data.frame(term_mat)
  
  df_aug <- df_lbl %>%
    mutate(
      effect_T       = term_df[["s(T)"]],
      effect_t_round = term_df[["s(t_round)"]],
      effect_RH      = term_df[["s(RH)"]],
      resid          = resid(gm),
      fitted         = fitted(gm),
      label          = lbl
    )
  
  # -------------------------------------------------------------------
  # Diagnostic plots: RAW residuals only
  # -------------------------------------------------------------------
  df_diag_raw <- df_aug %>%
    transmute(resid = resid, fitted, T, t_round, RH) %>%
    mutate(time_cal = origin_date + t_round)
  
  # residual vs fitted
  p1r <- ggplot(df_diag_raw, aes(x = fitted, y = resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, color = "red") +
    labs(x = "Fitted values (ppm)", y = "Residuals (ppm)", title = NULL) +
    base_theme_ieee
  
  # residual vs T
  p2r <- ggplot(df_diag_raw, aes(x = T, y = resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, color = "red") +
    labs(x = "T (\u00B0C)", y = "Residuals (ppm)", title = NULL) +
    base_theme_ieee
  
  # residual vs time (calendar axis, formatted as "YYYY-MM")
  p3r <- ggplot(df_diag_raw, aes(x = time_cal, y = resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, color = "red") +
    scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
    labs(x = "Time", y = "Residuals (ppm)", title = NULL) +
    base_theme_ieee +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
    )
  
  # residual vs RH
  p4r <- ggplot(df_diag_raw, aes(x = RH, y = resid)) +
    geom_point(alpha = 0.4) +
    geom_hline(yintercept = 0, color = "red") +
    labs(x = "RH (%)", y = "Residuals (ppm)", title = NULL) +
    base_theme_ieee
  
  # QQ-plot
  p5r <- ggplot(df_diag_raw, aes(sample = resid)) +
    stat_qq(alpha = 0.5) +
    stat_qq_line(color = "red") +
    labs(title = NULL, x = "Theoretical Q", y = "Sample Q") +
    base_theme_ieee
  
  # histogram
  p6r <- ggplot(df_diag_raw, aes(x = resid)) +
    geom_histogram(bins = 30, alpha = 0.8) +
    labs(title = NULL, x = "Residuals (ppm)", y = "Count") +
    base_theme_ieee
  
  # Save each plot (IEEE single-column size)
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_vs_fitted.png")),
    p1r, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_vs_T.png")),
    p2r, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_vs_time.png")),
    p3r, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_vs_RH.png")),
    p4r, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_qq.png")),
    p5r, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_hist.png")),
    p6r, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  
  # --------- ACF plot for irregular time index (t_round) ----------
  lags_unbiased <- acf_irregular(
    residuals = df_diag_raw$resid,
    t         = df_diag_raw$t_round,
    lag.max   = 14,
    unbiased  = TRUE
  )
  
  # make sure it's a plain data.frame
  lags_unbiased <- as.data.frame(lags_unbiased)
  
  # add lag-specific CI using number of pairs at each lag
  lags_unbiased$se       <- ifelse(lags_unbiased$pairs > 0,
                                   1.96 / sqrt(lags_unbiased$pairs),
                                   NA_real_)
  lags_unbiased$ci_upper <-  lags_unbiased$se
  lags_unbiased$ci_lower <- -lags_unbiased$se
  
  lags_no0 <- subset(lags_unbiased, lag != 0)
  
  ymax     <- max(abs(lags_no0$acf), na.rm = TRUE)
  ylim_val <- max(ymax, 0.05)
  
  p_acf <- ggplot() +
    geom_hline(yintercept = 0, color = "grey40") +
    geom_segment(
      data = lags_no0,
      aes(x = lag, xend = lag, y = acf, yend = 0)
    ) +
    geom_point(
      data = lags_no0,
      aes(x = lag, y = acf),
      size = 1
    ) +
    coord_cartesian(ylim = c(-ylim_val, ylim_val)) + 
    scale_x_continuous(breaks = 0:max(lags_no0$lag, na.rm = TRUE)) +
    labs(
      x = "Lag (days)",
      y = "Acf",
      title = NULL
    ) +
    base_theme_ieee
  
  ggsave(
    file.path(plot_dir, paste0(lbl, "_resid_acf.png")),
    p_acf, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  
  # -------------------------------------------------------------------
  # Smooth effects plots with confidence intervals
  # -------------------------------------------------------------------
  # typical values for other covariates
  typical_vals <- df_aug %>%
    summarise(
      t_round = median(t_round, na.rm = TRUE),
      RH      = median(RH,      na.rm = TRUE),
      T_med   = median(T,       na.rm = TRUE)
    )
  
  # grid for each covariate, holding others at typical values
  new_T <- tibble(
    T       = seq(min(df_aug$T, na.rm = TRUE),
                  max(df_aug$T, na.rm = TRUE),
                  length.out = 100),
    t_round = typical_vals$t_round,
    RH      = typical_vals$RH
  )
  
  new_t_round <- tibble(
    T       = typical_vals$T_med,
    t_round = seq(min(df_aug$t_round, na.rm = TRUE),
                  max(df_aug$t_round, na.rm = TRUE),
                  length.out = 100),
    RH      = typical_vals$RH
  )
  
  new_RH <- tibble(
    T       = typical_vals$T_med,
    t_round = typical_vals$t_round,
    RH      = seq(min(df_aug$RH, na.rm = TRUE),
                  max(df_aug$RH, na.rm = TRUE),
                  length.out = 100)
  )
  
  # predict terms + se for each grid
  pred_T       <- predict(gm, newdata = new_T,       type = "terms", se.fit = TRUE)
  pred_t_round <- predict(gm, newdata = new_t_round, type = "terms", se.fit = TRUE)
  pred_RH      <- predict(gm, newdata = new_RH,      type = "terms", se.fit = TRUE)
  
  effect_T  <- pred_T$fit[,"s(T)"]
  se_T      <- pred_T$se.fit[,"s(T)"]
  effect_tr <- pred_t_round$fit[,"s(t_round)"]
  se_tr     <- pred_t_round$se.fit[,"s(t_round)"]
  effect_RH <- pred_RH$fit[,"s(RH)"]
  se_RH     <- pred_RH$se.fit[,"s(RH)"]
  
  # assemble into data frames
  effects_T_df <- new_T %>%
    transmute(
      x      = T,
      effect = effect_T,
      lower  = effect_T - 1.96 * se_T,
      upper  = effect_T + 1.96 * se_T
    )
  
  # time grid as calendar date (instrumental drift)
  effects_tr_df <- new_t_round %>%
    transmute(
      t_round,
      time_cal = origin_date + t_round,
      effect   = effect_tr,
      lower    = effect_tr - 1.96 * se_tr,
      upper    = effect_tr + 1.96 * se_tr
    )
  
  effects_RH_df <- new_RH %>%
    transmute(
      x      = RH,
      effect = effect_RH,
      lower  = effect_RH - 1.96 * se_RH,
      upper  = effect_RH + 1.96 * se_RH
    )
  
  # ---- per-label plots ----
  common_effect_theme <- base_theme_ieee +
    theme(
      legend.position   = "bottom",
      legend.background = element_rect(fill = "transparent", colour = NA),
      legend.key        = element_blank()
    )
  
  # T effect
  p_T_effect <- ggplot(effects_T_df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI")) +
    geom_line(aes(y = effect, color = "Effect")) +
    geom_hline(aes(yintercept = 0, color = "y = 0 reference"), linetype = "dashed") +
    scale_fill_manual(name = NULL, values = c("95% CI" = "grey80")) +
    scale_color_manual(
      name   = NULL,
      values = c("Effect" = "black", "y = 0 reference" = "red")
    ) +
    labs(
      x = "T (\u00B0C)",
      y = "Effects (ppm)",
      title = NULL
    ) +
    common_effect_theme
  
  # RH effect
  p_RH_effect <- ggplot(effects_RH_df, aes(x = x)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI")) +
    geom_line(aes(y = effect, color = "Effect")) +
    geom_hline(aes(yintercept = 0, color = "y = 0 reference"), linetype = "dashed") +
    scale_fill_manual(name = NULL, values = c("95% CI" = "grey80")) +
    scale_color_manual(
      name   = NULL,
      values = c("Effect" = "black", "y = 0 reference" = "red")
    ) +
    labs(
      x = "RH (%)",
      y = "Effects (ppm)",
      title = NULL
    ) +
    common_effect_theme
  
  # Time (instrumental drift) effect, with calendar x-axis & denser ticks
  p_time_effect <- ggplot(effects_tr_df, aes(x = time_cal)) +
    geom_ribbon(aes(ymin = lower, ymax = upper, fill = "95% CI")) +
    geom_line(aes(y = effect, color = "Effect")) +
    geom_hline(aes(yintercept = 0, color = "y = 0 reference"), linetype = "dashed") +
    scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
    scale_fill_manual(name = NULL, values = c("95% CI" = "grey80")) +
    scale_color_manual(
      name   = NULL,
      values = c("Effect" = "black", "y = 0 reference" = "red")
    ) +
    labs(
      x = "Time",
      y = "Instrumental drift (ppm)",
      title = NULL
    ) +
    common_effect_theme +
    theme(
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
    )
  
  # save per-label effect plots
  ggsave(
    file.path(plot_dir, paste0(lbl, "_effect_T.png")),
    p_T_effect, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_effect_RH.png")),
    p_RH_effect, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  ggsave(
    file.path(plot_dir, paste0(lbl, "_effect_time.png")),
    p_time_effect, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  
  # -------------------------------------------------------------------
  # Fitted vs Observed plot (no title)
  # -------------------------------------------------------------------
  df_fit_obs <- df_aug %>%
    transmute(
      observed = err,
      fitted   = fitted
    )
  
  p_fit_obs <- ggplot(df_fit_obs, aes(x = fitted, y = observed)) +
    geom_point(alpha = 0.4) +
    geom_abline(slope = 1, intercept = 0, color = "red") +
    labs(
      x = "Fitted values (ppm)",
      y = "Observed values (ppm)",
      title = NULL
    ) +
    base_theme_ieee
  
  ggsave(
    file.path(plot_dir, paste0(lbl, "_fitted_vs_observed.png")),
    p_fit_obs, width = ieee_width, height = ieee_height, dpi = ieee_dpi
  )
  
  # -------------------------------------------------------------------
  # Model statistics + concurvity estimates
  # -------------------------------------------------------------------
  sg  <- summary(gm)
  
  conc_full <- concurvity(gm, full = TRUE)
  conc_est  <- conc_full["estimate", ]
  
  conc_T       <- as.numeric(conc_est["s(T)"])
  conc_t_round <- as.numeric(conc_est["s(t_round)"])
  conc_RH      <- as.numeric(conc_est["s(RH)"])
  
  stats <- tibble(
    label      = lbl,
    n_obs      = nrow(df_lbl),
    r_sq       = sg$r.sq,
    dev_expl   = sg$dev.expl,
    scale_est  = sg$sig2,
    gcv_ubre   = sg$gcv.ubre,
    edf_total  = sum(sg$edf),
    aic_gam    = AIC(gm),
    logLik_gam = as.numeric(logLik(gm)),
    mean_resid = mean(resid(gm)),
    sd_resid   = sd(resid(gm)),
    rmse       = sqrt(mean(resid(gm)^2)),
    
    conc_T       = conc_T,
    conc_t_round = conc_t_round,
    conc_RH      = conc_RH
  )
  
  all_augmented[[lbl]] <- df_aug
  model_stats[[lbl]]   <- stats
}

# -------------------------------------------------------------------
# Combine all per-label data
# -------------------------------------------------------------------
df_augmented_all <- bind_rows(all_augmented)
df_model_stats   <- bind_rows(model_stats)

# -------------------------------------------------------------------
# Build average smooths across labels for T, RH, and t_round
# -------------------------------------------------------------------

# Overall "typical" values for other covariates
typical_all <- df_augmented_all %>%
  summarise(
    t_round = median(t_round, na.rm = TRUE),
    RH      = median(RH,      na.rm = TRUE),
    T_med   = median(T,       na.rm = TRUE)
  )

# Common grids across all data
T_grid <- tibble(
  T       = seq(min(df_augmented_all$T, na.rm = TRUE),
                max(df_augmented_all$T, na.rm = TRUE),
                length.out = 200),
  t_round = typical_all$t_round,
  RH      = typical_all$RH
)

RH_grid <- tibble(
  T       = typical_all$T_med,
  t_round = typical_all$t_round,
  RH      = seq(min(df_augmented_all$RH, na.rm = TRUE),
                max(df_augmented_all$RH, na.rm = TRUE),
                length.out = 200)
)

t_round_grid <- tibble(
  T       = typical_all$T_med,
  t_round = seq(min(df_augmented_all$t_round, na.rm = TRUE),
                max(df_augmented_all$t_round, na.rm = TRUE),
                length.out = 200),
  RH      = typical_all$RH
)

# For each label-specific GAM, get partial effects on these grids
sT_mat <- map_dfc(models, function(m) {
  pred_terms <- predict(m, newdata = T_grid, type = "terms")
  pred_terms[,"s(T)"]
})

sRH_mat <- map_dfc(models, function(m) {
  pred_terms <- predict(m, newdata = RH_grid, type = "terms")
  pred_terms[,"s(RH)"]
})

sTR_mat <- map_dfc(models, function(m) {
  pred_terms <- predict(m, newdata = t_round_grid, type = "terms")
  pred_terms[,"s(t_round)"]
})

# Long data for all individual curves
effects_T_all <- purrr::map2_dfr(
  .x = models,
  .y = names(models),
  .f = function(m, lbl) {
    pred_terms <- predict(m, newdata = T_grid, type = "terms")
    tibble(
      label    = lbl,
      T        = T_grid$T,
      effect_T = pred_terms[, "s(T)"]
    )
  }
)

effects_RH_all <- purrr::map2_dfr(
  .x = models,
  .y = names(models),
  .f = function(m, lbl) {
    pred_terms <- predict(m, newdata = RH_grid, type = "terms")
    tibble(
      label     = lbl,
      RH        = RH_grid$RH,
      effect_RH = pred_terms[, "s(RH)"]
    )
  }
)

effects_time_all <- purrr::map2_dfr(
  .x = models,
  .y = names(models),
  .f = function(m, lbl) {
    pred_terms <- predict(m, newdata = t_round_grid, type = "terms")
    tibble(
      label     = lbl,
      t_round   = t_round_grid$t_round,
      time_cal  = origin_date + t_round_grid$t_round,
      effect_tr = pred_terms[, "s(t_round)"]
    )
  }
)

# Row-wise mean and SD across labels
avg_sT  <- rowMeans(as.matrix(sT_mat))
sd_sT   <- apply(as.matrix(sT_mat), 1, sd)

avg_sRH <- rowMeans(as.matrix(sRH_mat))
sd_sRH  <- apply(as.matrix(sRH_mat), 1, sd)

avg_sTR <- rowMeans(as.matrix(sTR_mat))
sd_sTR  <- apply(as.matrix(sTR_mat), 1, sd)

avg_T_smooth <- T_grid %>%
  transmute(
    T,
    avg_effect_T = avg_sT,
    sd_effect_T  = sd_sT
  )

avg_RH_smooth <- RH_grid %>%
  transmute(
    RH,
    avg_effect_RH = avg_sRH,
    sd_effect_RH  = sd_sRH
  )

avg_time_smooth <- t_round_grid %>%
  transmute(
    t_round,
    time_cal       = origin_date + t_round,
    avg_effect_tr  = avg_sTR,
    sd_effect_tr   = sd_sTR
  )

# -------------------------------------------------------------------
# Average effect plots
# -------------------------------------------------------------------

avg_common_theme <- base_theme_ieee +
  theme(
    legend.position   = "bottom",
    legend.background = element_rect(fill = "transparent", colour = NA),
    legend.key        = element_blank()
  )

# T
p_T_effects_avg <- ggplot() +
  geom_ribbon(
    data = avg_T_smooth,
    aes(
      x    = T,
      ymin = avg_effect_T - 2 * sd_effect_T,
      ymax = avg_effect_T + 2 * sd_effect_T,
      fill = "\u00b1 2 SD"
    )
  ) +
  geom_line(
    data = effects_T_all,
    aes(
      x     = T,
      y     = effect_T,
      group = label,
      color = "Individual"
    ),
    alpha     = 0.9,
    linewidth = 0.3
  ) +
  geom_line(
    data = avg_T_smooth,
    aes(x = T, y = avg_effect_T, color = "Average")
  ) +
  geom_hline(
    aes(yintercept = 0, color = "y = 0 reference"),
    linetype = "dashed"
  ) +
  scale_color_manual(
    name   = NULL,
    breaks = c("Average", "y = 0 reference", "Individual"),
    values = c(
      "Average"        = "black",
      "y = 0 reference"= "red",
      "Individual"     = "grey70"
    )
  ) +
  scale_fill_manual(
    name   = NULL,
    breaks = c("\u00b1 2 SD"),
    values = c("\u00b1 2 SD" = "grey80")
  ) +
  guides(
    color = guide_legend(order = 1),
    fill  = guide_legend(order = 2)
  ) +
  labs(
    x = "T (\u00B0C)",
    y = "Effects (ppm)",
    title = NULL
  ) +
  avg_common_theme

ggsave(
  file.path(plot_dir, "avg_effect_T.png"),
  p_T_effects_avg, width = ieee_width, height = ieee_height, dpi = ieee_dpi
)

# RH
p_RH_effects_avg <- ggplot() +
  geom_ribbon(
    data = avg_RH_smooth,
    aes(
      x    = RH,
      ymin = avg_effect_RH - 2 * sd_effect_RH,
      ymax = avg_effect_RH + 2 * sd_effect_RH,
      fill = "\u00b1 2 SD"
    )
  ) +
  geom_line(
    data = effects_RH_all,
    aes(
      x     = RH,
      y     = effect_RH,
      group = label,
      color = "Individual"
    ),
    alpha     = 0.9,
    linewidth = 0.3
  ) +
  geom_line(
    data = avg_RH_smooth,
    aes(x = RH, y = avg_effect_RH, color = "Average")
  ) +
  geom_hline(
    aes(yintercept = 0, color = "y = 0 reference"),
    linetype = "dashed"
  ) +
  scale_color_manual(
    name   = NULL,
    breaks = c("Average", "y = 0 reference", "Individual"),
    values = c(
      "Average"        = "black",
      "y = 0 reference"= "red",
      "Individual"     = "grey70"
    )
  ) +
  scale_fill_manual(
    name   = NULL,
    breaks = c("\u00b1 2 SD"),
    values = c("\u00b1 2 SD" = "grey80")
  ) +
  guides(
    color = guide_legend(order = 1),
    fill  = guide_legend(order = 2)
  ) +
  labs(
    x = "RH (%)",
    y = "Effects (ppm)",
    title = NULL
  ) +
  avg_common_theme

ggsave(
  file.path(plot_dir, "avg_effect_RH.png"),
  p_RH_effects_avg, width = ieee_width, height = ieee_height, dpi = ieee_dpi
)

# Time / instrumental drift
p_time_effects_avg <- ggplot() +
  geom_ribbon(
    data = avg_time_smooth,
    aes(
      x    = time_cal,
      ymin = avg_effect_tr - 2 * sd_effect_tr,
      ymax = avg_effect_tr + 2 * sd_effect_tr,
      fill = "\u00b1 2 SD"
    )
  ) +
  geom_line(
    data = effects_time_all,
    aes(
      x     = time_cal,
      y     = effect_tr,
      group = label,
      color = "Individual"
    ),
    alpha     = 0.9,
    linewidth = 0.3
  ) +
  geom_line(
    data = avg_time_smooth,
    aes(x = time_cal, y = avg_effect_tr, color = "Average")
  ) +
  geom_hline(
    aes(yintercept = 0, color = "y = 0 reference"),
    linetype = "dashed"
  ) +
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months") +
  scale_color_manual(
    name   = NULL,
    breaks = c("Average", "y = 0 reference", "Individual"),
    values = c(
      "Average"        = "black",
      "y = 0 reference"= "red",
      "Individual"     = "grey70"
    )
  ) +
  scale_fill_manual(
    name   = NULL,
    breaks = c("\u00b1 2 SD"),
    values = c("\u00b1 2 SD" = "grey80")
  ) +
  guides(
    color = guide_legend(order = 1),
    fill  = guide_legend(order = 2)
  ) +
  labs(
    x = "Time",
    y = "Instrumental drift (ppm)",
    title = NULL
  ) +
  avg_common_theme +
  theme(
    axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1)
  )

ggsave(
  file.path(plot_dir, "avg_effect_time.png"),
  p_time_effects_avg, width = ieee_width, height = ieee_height, dpi = ieee_dpi
)

# -------------------------------------------------------------------
# Map the average T and RH smooths back to each row in df_augmented_all
# -------------------------------------------------------------------
df_augmented_all <- df_augmented_all %>%
  mutate(
    effect_T_avg = approx(
      x    = avg_T_smooth$T,
      y    = avg_T_smooth$avg_effect_T,
      xout = T,
      rule = 2
    )$y,
    
    effect_RH_avg = approx(
      x    = avg_RH_smooth$RH,
      y    = avg_RH_smooth$avg_effect_RH,
      xout = RH,
      rule = 2
    )$y
  )

# Save data and models
write_csv(df_augmented_all, "./saved_data/R/gam_augmented_all_labels.csv")
write_csv(df_model_stats,   "./saved_data/R/gam_model_stats_by_label.csv")
saveRDS(models, file = "./saved_data/R/gam_models_list.rds")

# quick preview
print(df_model_stats)

# -------------------------------------------------------------------
# Create simple CSV with R2, RMSE, normalized RMSE (RMSE / SD)
# -------------------------------------------------------------------

# Load saved tables (if running fresh)
df_stats <- read_csv("./saved_data/R/gam_model_stats_by_label.csv",
                     show_col_types = FALSE)

df_aug_all <- read_csv("./saved_data/R/gam_augmented_all_labels.csv",
                       show_col_types = FALSE)

# Compute SD of observed err per label
sd_by_label <- df_aug_all %>%
  group_by(label) %>%
  summarise(
    err_sd = sd(err, na.rm = TRUE),
    .groups = "drop"
  )

# Combine with model stats
simple_metrics <- df_stats %>%
  select(label, r_sq, rmse) %>%
  left_join(sd_by_label, by = "label") %>%
  mutate(
    R1       = r_sq,
    RMSE     = rmse,
    nRMSE_sd = if_else(err_sd > 0, RMSE / err_sd, NA_real_)
  ) %>%
  select(label, R1, RMSE, nRMSE_sd) %>%
  # round numeric columns to 3 decimals
  mutate(across(c(R1, RMSE, nRMSE_sd), ~ round(.x, 3)))

# Save CSV
write_csv(simple_metrics,
          "./saved_data/R/gam_simple_metrics_by_label.csv")

#### 
concurvity_metrics <- df_stats %>%
  select(label, conc_T, conc_t_round, conc_RH) %>%
  mutate(
    across(
      c(conc_T, conc_t_round, conc_RH),
      ~ round(.x, 3)
    )
  )

write_csv(
  concurvity_metrics,
  "./saved_data/R/gam_concurvity_metrics_by_label.csv"
)
