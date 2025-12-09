rm(list = ls())
library(mgcv)
library(dplyr)
library(readr)
library(ggplot2)
library(patchwork)

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

data <- read_csv("./saved_data/R/gam_data_outlier_removed.csv", show_col_types = FALSE)

# label_list <- c("1187", "1215", "1029", "1045", "1312", "1288", "1022",
#                 "1097", "1098", "1094", "1103", "1104", "1055")
label_list <- c("1097") # for testing

## helper to compute stats ----------------------------------------------
model_stats <- function(mod, data, response = "err") {
  y    <- data[[response]]
  yhat <- predict(mod, newdata = data)
  
  # RMSE
  rmse <- sqrt(mean((y - yhat)^2, na.rm = TRUE))
  
  # R² (from predictions)
  ss_res <- sum((y - yhat)^2, na.rm = TRUE)
  ss_tot <- sum((y - mean(y, na.rm = TRUE))^2, na.rm = TRUE)
  r2 <- if (ss_tot == 0) NA_real_ else 1 - ss_res / ss_tot
  
  # Standard deviation of response
  sd_y <- sd(y, na.rm = TRUE)
  rmse_std <- if (sd_y == 0 || is.na(sd_y)) NA_real_ else rmse / sd_y
  
  tibble(
    R2        = r2,
    RMSE      = rmse,
    RMSE_std  = rmse_std,
    AIC       = AIC(mod)
  )
}

## object to store results ----------------------------------------------
results <- tibble(
  label    = character(),
  model    = character(),
  R2       = numeric(),
  RMSE     = numeric(),
  RMSE_std = numeric(),
  AIC      = numeric()
)

## object to store fitted models ----------------------------------------
all_models <- list()

## loop over labels ------------------------------------------------------
for (lbl in label_list) {
  message("Processing label: ", lbl)
  start_time <- Sys.time()
  
  # Subset data for this label
  df_lbl <- dplyr::filter(data, label == lbl)
  
  # -------------------------------------------------------------------
  # Model 1: baseline (no interaction)
  # -------------------------------------------------------------------
  gm1 <- bam(
    err ~ s(T, k = 10, bs = "cr") +
      s(t_round, k = 200, bs = "cr", pc = 0) +
      s(RH, k = 10, bs = "cr"),
    data = df_lbl, method = "fREML", discrete = TRUE
  )
  message("Model 1 fitting finished.")
  
  stats_gm1 <- model_stats(gm1, df_lbl) |>
    mutate(label = lbl, model = "baseline_T + t_round + RH")
  
  # -------------------------------------------------------------------
  # Model 2: interaction of T and RH + t_round effect
  # -------------------------------------------------------------------
  gm2 <- bam(
    err ~ 
      s(T, k = 10, bs = "cr") +
      s(t_round, k = 200, bs = "cr", pc = 0) +
      s(RH, k = 10, bs = "cr") +
      ti(T, RH, k = c(10, 10), bs = c("cr", "cr")),
    data = df_lbl, method = "fREML", discrete = TRUE
  )
  message("Model 2 fitting finished.")
  
  stats_gm2 <- model_stats(gm2, df_lbl) |>
    mutate(label = lbl, model = "s(T) + s(RH) + ti(T,RH) + s(t_round)")
  
  # -------------------------------------------------------------------
  # Model 3: all pairwise interactions (T–RH, T–t_round, RH–t_round)
  # -------------------------------------------------------------------
  gm3 <- bam(
    err ~ 
      s(T, k = 10, bs = "cr") +
      s(t_round, k = 200, bs = "cr", pc = 0) +
      s(RH, k = 10, bs = "cr") +
      ti(T, RH,      k = c(10, 10),  bs = c("cr", "cr")) +
      ti(T, t_round, k = c(10, 200), bs = c("cr", "cr")) +
      ti(RH, t_round,k = c(10, 200), bs = c("cr", "cr")),
    data = df_lbl, method = "fREML", discrete = TRUE
  )
  message("Model 3 fitting finished.")
  
  stats_gm3 <- model_stats(gm3, df_lbl) |>
    mutate(label = lbl, model = "s(T)+s(RH)+s(t_round)+all_pairwise_ti")
  
  # -------------------------------------------------------------------
  # Model 4: full model with three-way interaction ti(T, RH, t_round)
  # -------------------------------------------------------------------
  gm4 <- bam(
    err ~ 
      s(T, k = 10, bs = "cr") +
      s(t_round, k = 200, bs = "cr", pc = 0) +
      s(RH, k = 10, bs = "cr") +
      ti(T, RH,      k = c(10, 10),   bs = c("cr", "cr")) +
      ti(T, t_round, k = c(10, 200),  bs = c("cr", "cr")) +
      ti(RH, t_round,k = c(10, 200),  bs = c("cr", "cr")) +
      ti(T, RH, t_round, k = c(5, 5, 50), bs = c("cr", "cr", "cr")),
    data = df_lbl, method = "fREML", discrete = TRUE
  )
  message("Model 4 (with three-way interaction) fitting finished.")
  
  stats_gm4 <- model_stats(gm4, df_lbl) |>
    mutate(label = lbl, model = "s(T)+s(RH)+s(t_round)+all_pairwise_ti+ti(T,RH,t_round)")
  
  # bind all stats for this label
  stats_lbl <- bind_rows(
    stats_gm1,
    stats_gm2,
    stats_gm3,
    stats_gm4
  )
  
  # print stats for this label
  message("Statistics for label ", lbl, ":")
  print(stats_lbl)
  
  # compute and print elapsed time
  elapsed <- Sys.time() - start_time
  elapsed_min <- as.numeric(elapsed, units = "secs") / 60
  message(sprintf("Time spent on label %s: %.2f minutes", lbl, elapsed_min))
  
  # add to overall results
  results <- bind_rows(results, stats_lbl)
  
  # save fitted models for this label into all_models
  all_models[[lbl]] <- list(
    gm1 = gm1,
    gm2 = gm2,
    gm3 = gm3,
    gm4 = gm4
  )
}

results <- results %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))

## save results and models to file ---------------------------------------
write_csv(results, "./saved_data/R/gam_with_interaction_comparison_bam.csv")
saveRDS(all_models, "./saved_data/R/gam_with_interaction_models_bam.rds")



