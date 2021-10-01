########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Script for preparing calibration data
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# People, Place, and Health Collective, Department of Epidemiology, Brown University
#
# Created: Dec 30, 2020
# Last update: June 23, 2021
#
# Purpose: Prepares all data required for calibration script
#
##################################################################################
###################        Cost effectiveness analysis      ######################
##################################################################################


# Load required packages
library(dplyr)
library(abind)
library(FME)
library(openxlsx)

source("data_input.R")

# INPUT PARAMETERS
sw.EMS.ODloc <- "overall" # Please choose from "overall" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "overall"

## Model parameter updates for calibration process ##
params <- data_input()
CalibPar <- read.xlsx(params$workbook, sheet = "CalibPar")
parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
row.names(parRange) <- CalibPar$par
set.seed(5112021)
calib.par <- Latinhyper(parRange, sample.size)
calib.par <- data.frame(calib.par)
saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))
# }
toc()

cal_param_names <- names(calib.par)

calib.parameters <- list()

tic("outer loop")
for (bb in 1:(sample.size / batch.size)) {
  batch.bgn <- (bb - 1) * batch.size + 1
  batch.end <- bb * batch.size
  ii <- 1
  tic("inner loop")
  for (cc in batch.bgn:batch.end) {
    for (pp in 1:length(cal_param_names)) {
      params[[cal_param_names[pp]]] <- calib.par[cc, pp]
    }
    # Overdose probability matrix (per month)
    overdose_probs <- matrix(0, nrow = 4, ncol = 2)
    rownames(overdose_probs) <- c("rx", "il_lr", "il_hr", "NODU")
    colnames(overdose_probs) <- c("first", "subs")
    overdose_probs["rx", "subs"] <- params$od_rx_sub
    overdose_probs["il_lr", "subs"] <- params$od_il_lr_sub
    overdose_probs["il_hr", "subs"] <- params$od_il_lr_sub * params$multi_hr
    overdose_probs["NODU", "subs"] <- params$od_nodu_sub
    overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi_sub
    params$overdose_probs <- overdose_probs

    # Baseline mortality excluding overdose (per month)
    mortality_probs <- matrix(0, nrow = 2, ncol = length(mor.gp))
    rownames(mortality_probs) <- c("bg", "drug")
    colnames(mortality_probs) <- mor.gp
    mortality_probs["bg", ] <- params$mortality_base
    mortality_probs["drug", ] <- params$mortality_drug
    params$mortality_probs <- mortality_probs
    params$od_911_pub <- params$od_911_priv * params$od_911_pub_mul

    calib.parameters[[ii]] <- params
    ii <- ii + 1
  }
  toc()
  saveRDS(calib.parameters, file = paste0("Inputs/CalibrationSampleDatax", bb, ".rds"))
}
toc()


## Group into batches
calib.par$batch_number <- (as.integer(rownames(calib.par)) - 1) / (batch.size) + 1
calib.par$batch_number <- trunc(calib.par$batch_number)

prep_data <- function(cc) {
  for (pp in 1:length(cal_param_names)) {
    params[[cal_param_names[pp]]] <- calib.par[cc, pp]
  }
  ii <- cc - batch.bgn + 1
  # Overdose probability matrix (per month)
  overdose_probs <- matrix(0, nrow = 4, ncol = 2)
  rownames(overdose_probs) <- c("rx", "il_lr", "il_hr", "NODU")
  colnames(overdose_probs) <- c("first", "subs")
  overdose_probs["rx", "subs"] <- params$od_rx_sub
  overdose_probs["il_lr", "subs"] <- params$od_il_lr_sub
  overdose_probs["il_hr", "subs"] <- params$od_il_lr_sub * params$multi_hr
  overdose_probs["NODU", "subs"] <- params$od_nodu_sub
  overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi_sub
  params$overdose_probs <- overdose_probs

  # Baseline mortality excluding overdose (per month)
  mortality_probs <- matrix(0, nrow = 2, ncol = length(mor.gp))
  rownames(mortality_probs) <- c("bg", "drug")
  colnames(mortality_probs) <- mor.gp
  mortality_probs["bg", ] <- params$mortality_base
  mortality_probs["drug", ] <- params$mortality_drug
  params$mortality_probs <- mortality_probs
  params$od_911_pub <- params$od_911_priv * params$od_911_pub_mul

  calib.parameters[[ii]] <- params
}
