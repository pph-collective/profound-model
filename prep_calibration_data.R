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
library(tictoc)

source("data_input.R")
# # INPUT PARAMETERS
# sw.EMS.ODloc <- "overall" # Please choose from "overall" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "overall"

## Model parameter updates for calibration process ##

CalibPar <- read.xlsx(WB, sheet = "CalibPar")
parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
row.names(parRange) <- CalibPar$par
set.seed(5112021)
calib.par <- Latinhyper(parRange, sample.size)
calib.par <- data.frame(calib.par)
saveRDS(calib.par, paste0("calibration/prep_calibration_data/Calib_par_table.rds"))

## Load precalibrated decision tree data 
precalib.par <- readRDS(file = "Inputs/Precalib_dtree.rds")
OD_pub <- precalib.par[,"OD_pub"]
rr_OD_wit_priv <- precalib.par[,"rr_OD_wit_priv"]
rr_OD_911_priv <- precalib.par[,"rr_OD_911_priv"]
  
# Initialize calibration parameters
cal_param_names <- names(calib.par)
calib.parameters <- list()

tic("outer loop")
for (bb in 1:(sample.size / batch.size)) {
  batch.bgn <- (bb - 1) * batch.size + 1
  batch.end <- bb * batch.size
  ii <- 1
  for (cc in batch.bgn:batch.end) {
    for (pp in 1:length(cal_param_names)) {
      params[[cal_param_names[pp]]] <- calib.par[cc, pp]
    }
    # Overdose probability matrix (per month)
    overdose_probs <- matrix(0, nrow = 4, ncol = 2)
    rownames(overdose_probs) <- c("preb", "il.lr", "il.hr", "NODU")
    colnames(overdose_probs) <- c("first", "subs")
    overdose_probs["preb", "subs"] <- params$od.preb.sub
    overdose_probs["il.lr", "subs"] <- params$od.il.lr.sub
    overdose_probs["il.hr", "subs"] <- params$od.il.lr.sub * params$multi.hr
    overdose_probs["NODU", "subs"] <- params$od.NODU.sub
    overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi.sub
    params$overdose_probs <- overdose_probs

    # Baseline mortality excluding overdose (per month)
    mortality_probs <- matrix(0, nrow = 2, ncol = length(mor.gp))
    rownames(mortality_probs) <- c("bg", "drug")
    colnames(mortality_probs) <- mor.gp
    mortality_probs["bg", ] <- params$mor.bg
    mortality_probs["drug", ] <- params$mor.drug
    params$mortality_probs <- mortality_probs
    params$OD_loc_pub  <- OD_pub[cc]
    params$OD_wit_priv <- params$OD_wit_pub * rr_OD_wit_priv[cc] 
    params$OD_911_priv <- params$OD_911_pub * rr_OD_911_priv[cc]

    calib.parameters[[ii]] <- params
    ii <- ii + 1
  }
  saveRDS(calib.parameters, file = paste0("calibration/prep_calibration_data/CalibrationSampleData", bb, ".rds"))
}
toc()


# ## Group into batches
# calib.par$batch_number <- (as.integer(rownames(calib.par)) - 1) / (batch.size) + 1
# calib.par$batch_number <- trunc(calib.par$batch_number)
# 
# prep_data <- function(cc) {
#   for (pp in 1:length(cal_param_names)) {
#     params[[cal_param_names[pp]]] <- calib.par[cc, pp]
#   }
#   ii <- cc - batch.bgn + 1
#   # Overdose probability matrix (per month)
#   overdose_probs <- matrix(0, nrow = 4, ncol = 2)
#   rownames(overdose_probs) <- c("preb", "il.lr", "il.hr", "NODU")
#   colnames(overdose_probs) <- c("first", "subs")
#   overdose_probs["preb", "subs"] <- params$od.preb.sub
#   overdose_probs["il.lr", "subs"] <- params$od.il.lr.sub
#   overdose_probs["il.hr", "subs"] <- params$od.il.lr.sub * params$multi.hr
#   overdose_probs["NODU", "subs"] <- params$od.NODU.sub
#   overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi.sub
#   params$overdose_probs <- overdose_probs
# 
#   # Baseline mortality excluding overdose (per month)
#   mortality_probs <- matrix(0, nrow = 2, ncol = length(mor.gp))
#   rownames(mortality_probs) <- c("bg", "drug")
#   colnames(mortality_probs) <- mor.gp
#   mortality_probs["bg", ] <- params$mor.bg
#   mortality_probs["drug", ] <- params$mor.drug
#   params$mortality_probs <- mortality_probs
#   params$OD_911_pub <- params$OD_911_priv * params$OD_911_pub_mul
# 
#   calib.parameters[[ii]] <- params
# }
