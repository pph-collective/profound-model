################################################################################################
#########################       Calibration data preparation     ###############################
################################################################################################

################################################################################################
# This is to prepare all the data required by the calibration procedure
# Load random samples for model parameters (based on latin hypercube sampling) and combine other parameters (uncalibrated) and save as parameter lists

#Load required packages
library(dplyr)
library(abind)
#INPUT PARAMETERS
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

source("Profound-DataInput.R")

## Model parameter updates for calibration process ##
if(file.exists(paste0("Inputs/Calib_par_table.rds"))){
  calib.par  <- readRDS(paste0("Inputs/Calib_par_table.rds"))
} else if (!file.exists(paste0("Inputs/Calib_par_table.rds"))){
  library(FME)
  CalibPar <- read.xlsx(WB, sheet="CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par
  set.seed(5112021)
  calib.par <- Latinhyper(parRange, sample.size)
  calib.par <- data.frame(calib.par)
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))
}

nm.calp <- names(calib.par)

calib.parameters <- list()

for (bb in 1:(sample.size/batch.size)){
  batch.bgn <- (bb-1)*batch.size + 1
  batch.end <- bb*batch.size
  ii <- 1
  for (cc in batch.bgn:batch.end){
    for (pp in 1:length(nm.calp)){
      vparameters[[nm.calp[pp]]] <- calib.par[cc, pp]
    }
    # Overdose probability matrix (per month)
    od.matrix             <- matrix(0, nrow = 4, ncol = 2)
    rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
    colnames(od.matrix)   <- c("first", "subs")
    od.matrix["preb", "subs"]   <- vparameters$od.preb.sub
    od.matrix["il.lr", "subs"]  <- vparameters$od.il.lr.sub
    od.matrix["il.hr", "subs"]  <- vparameters$od.il.lr.sub * vparameters$multi.hr
    od.matrix["NODU", "subs"]   <- vparameters$od.NODU.sub
    od.matrix[ , "first"]       <- od.matrix[ , "subs"] / vparameters$multi.sub
    vparameters$od.matrix       <- od.matrix
    
    # Baseline mortality excluding overdose (per month)
    mor.matrix                  <- matrix(0, nrow = 2, ncol = length(mor.gp))
    rownames(mor.matrix)        <- c("bg", "drug")
    colnames(mor.matrix)        <- mor.gp
    mor.matrix["bg", ]          <- vparameters$mor.bg
    mor.matrix["drug", ]        <- vparameters$mor.drug
    vparameters$mor.matrix      <- mor.matrix
    vparameters$OD_911_pub      <- vparameters$OD_911_priv * vparameters$OD_911_pub_mul
    
    calib.parameters[[ii]] <- vparameters
    ii <- ii+1
  }
  saveRDS(calib.parameters, file = paste0("Inputs/CalibrationSampleData", bb, ".rds"))
}

