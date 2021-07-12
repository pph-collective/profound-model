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
##################################################################################
###################        Cost effectiveness analysis      ######################
##################################################################################


#Load required packages
library(dplyr)
library(abind)
library(FME)
library(openxlsx)

source("data_input.R")
#INPUT PARAMETERS
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"
sample.size <- 10
batch.size <- 5


## Model parameter updates for calibration process ##
tic("parameter updates")
# if(file.exists(paste0("Inputs/Calib_par_table.rds"))){
#   calib.par  <- readRDS(paste0("Inputs/Calib_par_table.rds"))
# } else {
  CalibPar <- read.xlsx(WB, sheet="CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par
  set.seed(5112021)
  calib.par <- Latinhyper(parRange, sample.size)
  calib.par <- data.frame(calib.par)
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))
# }
toc()

nm.calp <- names(calib.par)

calib.parameters <- list()

tic("outer loop")
for (bb in 1:(sample.size/batch.size)){
  batch.bgn <- (bb-1)*batch.size + 1
  batch.end <- bb*batch.size
  ii <- 1
  tic("inner loop")
  for (cc in batch.bgn:batch.end){
    stopifnot(ii + batch.bgn - 1 == cc)
    stopifnot(ii == cc - batch.bgn + 1)
    for (pp in 1:length(nm.calp)){
      params[[nm.calp[pp]]] <- calib.par[cc, pp]
    }
    # Overdose probability matrix (per month)
    od.matrix             <- matrix(0, nrow = 4, ncol = 2)
    rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
    colnames(od.matrix)   <- c("first", "subs")
    od.matrix["preb", "subs"]   <- params$od.preb.sub
    od.matrix["il.lr", "subs"]  <- params$od.il.lr.sub
    od.matrix["il.hr", "subs"]  <- params$od.il.lr.sub * params$multi.hr
    od.matrix["NODU", "subs"]   <- params$od.NODU.sub
    od.matrix[ , "first"]       <- od.matrix[ , "subs"] / params$multi.sub
    params$od.matrix       <- od.matrix
    
    # Baseline mortality excluding overdose (per month)
    mor.matrix                  <- matrix(0, nrow = 2, ncol = length(mor.gp))
    rownames(mor.matrix)        <- c("bg", "drug")
    colnames(mor.matrix)        <- mor.gp
    mor.matrix["bg", ]          <- params$mor.bg
    mor.matrix["drug", ]        <- params$mor.drug
    params$mor.matrix      <- mor.matrix
    params$OD_911_pub      <- params$OD_911_priv * params$OD_911_pub_mul
    
    calib.parameters[[ii]] <- params
    ii <- ii+1
  }
  toc()
  saveRDS(calib.parameters, file = paste0("Inputs/CalibrationSampleDatax", bb, ".rds"))
}
toc()


## Group into batches
calib.par$batch_number <- (as.integer(rownames(calib.par)) - 1) / (batch.size) + 1
calib.par$batch_number <- trunc(calib.par$batch_number)

prep_data <- function(cc) {
  for (pp in 1:length(nm.calp)){
      params[[nm.calp[pp]]] <- calib.par[cc, pp]
    }
    ii <- cc - batch.bgn + 1
    # Overdose probability matrix (per month)
    od.matrix             <- matrix(0, nrow = 4, ncol = 2)
    rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
    colnames(od.matrix)   <- c("first", "subs")
    od.matrix["preb", "subs"]   <- params$od.preb.sub
    od.matrix["il.lr", "subs"]  <- params$od.il.lr.sub
    od.matrix["il.hr", "subs"]  <- params$od.il.lr.sub * params$multi.hr
    od.matrix["NODU", "subs"]   <- params$od.NODU.sub
    od.matrix[ , "first"]       <- od.matrix[ , "subs"] / params$multi.sub
    params$od.matrix       <- od.matrix
    
    # Baseline mortality excluding overdose (per month)
    mor.matrix                  <- matrix(0, nrow = 2, ncol = length(mor.gp))
    rownames(mor.matrix)        <- c("bg", "drug")
    colnames(mor.matrix)        <- mor.gp
    mor.matrix["bg", ]          <- params$mor.bg
    mor.matrix["drug", ]        <- params$mor.drug
    params$mor.matrix      <- mor.matrix
    params$OD_911_pub      <- params$OD_911_priv * params$OD_911_pub_mul
    
    calib.parameters[[ii]] <- params
  
}
