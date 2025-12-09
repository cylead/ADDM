rm(list = ls())
library(readr)
library(tidyverse) 
# ---- Setup ----
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("functions.R")
data <- read_csv("./saved_data/R/gam_data_for_R.csv")

label_thresholds <- tribble(
  ~label, ~thresh_RH, ~thresh_err,
  "1187", 90, 4000,
  "1215", 90, 1500,
  "1029", 80, 2700,
  "1045", 80, 4500,
  "1312", 90, 2500,
  "1288", 90, 1500,
  "1022", 90, 500,
  "1097", 90, 600,
  "1098", 90, 600,
  "1094", 90, 2700,
  "1103", 75, 160,
  "1104", 80, 5700,
  "1055", 75, 90
) 

filtered_data <- tibble()  # empty container

for (i in seq_len(nrow(label_thresholds))) {
  lbl     <- label_thresholds$label[i]
  thr_rh  <- label_thresholds$thresh_RH[i]
  thr_err <- label_thresholds$thresh_err[i]
  
  df <- data %>%
    filter(label == lbl, RH < thr_rh, err < thr_err)
  
  # ------------------------------------------------------------------------------
  # Round timestamp to days and compute t_round (in "round_unit" units)
  # (Assumes helper functions get_rounder() and get_unit_factor() exist.)
  # ------------------------------------------------------------------------------
  round_unit <- "day"
  round_fun  <- "floor"
  
  rounder <- get_rounder(round_fun)
  
  df <- df %>%
    mutate(ts_round = rounder(timestamp, unit = round_unit))
  
  unit_factor <- get_unit_factor(round_unit)  # seconds per "round_unit"
  df <- df %>%
    mutate(t_round = as.numeric(difftime(ts_round, min(ts_round), units = "secs")) / unit_factor)
  
  
  filtered_data <- bind_rows(filtered_data, df)
}

# ---- Save result ----
write_csv(filtered_data, "./saved_data/R/gam_data_outlier_removed.csv")

