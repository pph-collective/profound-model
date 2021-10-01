#!/usr/bin/env Rscript

#' Analyze calibration and select best-fitting runs
#'
#' @description
#' imports model runs and finds the runs that fit empirical data best, for use in
#' future analyses
#'
#' @param TODO need to make function
#'
#' @returns
#' writes parameters for calibrated runs to file
#'

rm(list = ls())
# Load required packages
library(dplyr)
library(openxlsx)
library(abind)
library(ggplot2)
library(xlsx)

# Load scripts
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")
source("parallel.R")
source("data_input.R")


calibration_results <- NULL
# read in and bind calibration results
for (batch.ind in 1:10) {
  temp_results <- readRDS(paste0("calibration/CalibrationOutputs", batch.ind, ".rds"))
  calibration_results <- rbind(calibration_results, temp_results)
}
rm(temp_results)

# add the calibration parameters to the calibration results
calibration_params <- readRDS(paste0("calibration/Calib_par_table.rds"))
calibration_results <- cbind(calibration_results, calibration_params)

# read in workbook of targets
workbook <- loadWorkbook("Inputs/MasterTable.xlsx")
Target <- read.xlsx(workbook, sheet = "Target")
tar.data <- Target$pe


# Calculate goodness of fit (gof)
for (ss in 1:nrow(calibration_results)) {
  prediction <- calibration_results[ss, c(
    "od.death16", "od.death17", "od.death18", "od.death19",
    "fx.death16", "fx.death17", "fx.death18", "fx.death19",
    "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19"
  )]
  gof <- 0
  
  # TODO needs to be less hardcoded at some point
  for (j in 1:length(tar.data)) {
    if (j %in% c(5:8)) { # 5:8 is proportion that involves fx
      gof <- gof + (abs(prediction[j] * 100 - tar.data[j] * 100) / (tar.data[j] * 100)) / length(tar.data)
    } else {
      gof <- gof + (abs(prediction[j] - tar.data[j]) / tar.data[j]) / length(tar.data)
    }
  }
  calibration_results[ss, "gof"] <- gof
}

# sort by goodness of fit and select top 100 fits
sorted.mx <- calibration_results[order(calibration_results[, "gof"], decreasing = F), ]
cal_sample <- 100
calibration_results <- sorted.mx[1:cal_sample, ]
write.xlsx(calibration_results,
  file = paste0("calibration/Calibrated_results.xlsx"),
  col.names = T, row.names = F
)


## save calibrated results as parameter lists (prepare for main analysis)##
# INPUT PARAMETERS
sw.EMS.ODloc <- "overall" # Please choose from "overgit all" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "overall"
cal_param_names <- names(calibration_params)
calibrated_parameters <- list()


# for each selected calibration run, determine the run's parameters
for (cc in 1:nrow(calibration_results)) {
  for (pp in 1:length(cal_param_names)) {
    params[[cal_param_names[pp]]] <- calibration_results[cc, cal_param_names[pp]]
  }
  # Overdose probability parameter matrix (per month)
  overdose_probs <- matrix(0, nrow = 4, ncol = 2)
  rownames(overdose_probs) <- c("rx", "il_lr", "il_hr", "NODU")
  colnames(overdose_probs) <- c("first", "subs")
  overdose_probs["rx", "subs"] <- params$od_rx_sub
  overdose_probs["il_lr", "subs"] <- params$od_il_lr_sub
  overdose_probs["il_hr", "subs"] <- params$od_il_lr_sub * params$multi_hr
  overdose_probs["NODU", "subs"] <- params$od_nodu_sub
  overdose_probs[, "first"] <- overdose_probs[, "subs"] / params$multi_sub
  params$overdose_probs <- overdose_probs

  # Baseline mortality parameters, excluding overdose (per month)
  mortality_probs <- matrix(0, nrow = 2, ncol = length(mor.gp))
  rownames(mortality_probs) <- c("bg", "drug")
  colnames(mortality_probs) <- mor.gp
  mortality_probs["bg", ] <- params$mortality_base
  mortality_probs["drug", ] <- params$mortality_drug
  params$mortality_probs <- mortality_probs
  params$od_911_pub <- params$od_911_priv * params$od_911_pub_mul

  calibrated_parameters[[cc]] <- params
}

# alibrated.seed <- calibration_results[,"seed"]
saveRDS(calibrated_parameters, file = paste0("calibration/CalibratedData.rds"))
saveRDS(calibration_results[, "seed"], file = paste0("calibration/CalibratedSeed.rds"))



## Plot ##
par(mfrow = c(1, 3))
# REVIEWED md = median, change to mean
mean_oddeath <- apply(calibration_results[, c("od.death16", "od.death17", "od.death18", "od.death19")], 2, mean)
ymax <- max(calibration_results[, c("od.death16", "od.death17", "od.death18", "od.death19")])
plot(
  x = 2016:2019, tar.data[1:4], col = "black", pch = 18, xlab = "Year", ylab = "Number of opioid overdose deaths", cex = 1.2, cex.axis = 1.2, cex.lab = 1.3,
  xaxt = "n", ylim = c(0, 500), frame.plot = FALSE
)
# title("A", adj = 0, line =-0.5)
for (i in 1:cal_sample) {
  lines(x = 2016:2019, y = calibration_results[i, c("od.death16", "od.death17", "od.death18", "od.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd = 3)
}
lines(x = 2016:2019, y = mean_oddeath, col = "red", lwd = 4)
points(x = 2016:2019, tar.data[1:4], col = "black", pch = 16, cex = 1.2, cex.axis = 0.95)
axis(1, at = 2016:2019, pos = 0, lwd.ticks = 0, cex.axis = 1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h = 0)
legend("top",
  legend = c("Target", "Model: mean", "Model: simulation"),
  col = c("black", "red", adjustcolor("indianred1", alpha = 0.2)),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = 3,
  bty = "n",
  pt.cex = 1.2,
  cex = 1.1,
  text.col = "black",
  horiz = T
)
# legend(x="bottomright",
#        col=c("dodgerblue","firebrick2", "lightgrey"),
#        lwd=c(1.2, 1, 0.5), lty = c(1, NA, NA), pch=c(16, 16,16), pt.cex = c(1,1.1,1), cex = 0.8, bty = "n")


mean_fxdeath <- apply(calibration_results[, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], 2, median)
# REVIEWED: change from hardcoded ylim; mean instead of median for everything
plot(
  x = 2016:2019, tar.data[5:8], col = "black", pch = 18, xlab = "Year", ylab = "Proportion of OOD deaths with fentanyl present", cex = 1.2, cex.axis = 1.2, cex.lab = 1.3,
  xaxt = "n", yaxt = "n", ylim = c(0, 1), frame.plot = FALSE
)
# title("A", adj = 0, line =-0.5)
for (i in 1:cal_sample) {
  lines(x = 2016:2019, y = calibration_results[i, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd = 2)
}
lines(x = 2016:2019, y = mean_fxdeath, col = "red", lwd = 3)
points(x = 2016:2019, tar.data[5:8], col = "black", pch = 16, cex = 1.2, cex.axis = 0.95)
axis(1, at = 2016:2019, pos = 0, lwd.ticks = 0, cex.axis = 1.2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20%", "40%", "60%", "80%", "100%"), cex.axis = 1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h = 0)


md.edvisits <- apply(calibration_results[, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], 2, median)
plot(
  x = 2016:2019, tar.data[9:12], col = "black", pch = 18, xlab = "Year", ylab = "Number of ED visists for opioid overdose", cex = 1.2, cex.axis = 1.2, cex.lab = 1.3,
  xaxt = "n", ylim = c(0, 2500), frame.plot = FALSE
)

for (i in 1:cal_sample) {
  lines(x = 2016:2019, y = calibration_results[i, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], col = adjustcolor("indianred1", alpha = 0.2), lwd = 2)
}
lines(x = 2016:2019, y = md.edvisits, col = "red", lwd = 3)
points(x = 2016:2019, tar.data[9:12], col = "black", pch = 16, cex = 1.2, cex.axis = 1.2)
axis(1, at = 2016:2019, pos = 0, lwd.ticks = 0, cex.axis = 1.2)
abline(h = 0)

par(mfrow = c(1, 1))
legend("bottom",
  legend = c("Target", "Model: mean", "Model: simulation"),
  col = c("black", "red", adjustcolor("indianred1", alpha = 0.2)),
  pch = c(16, NA, NA),
  lty = c(NA, 1, 1),
  lwd = 3,
  bty = "n",
  pt.cex = 2,
  cex = 1.1,
  text.col = "black",
  horiz = T
)


cal_param_names <- names(calibration_params)


calib.post <- calibration_results[, (dim(calibration_results)[2] - length(cal_param_names) + 1):dim(calibration_results)[2]]
# REVIEWED par -> param
par <- rep(colnames(calib.post), 2)
# REVIEWED pe = point estimate, lend = lower end, uend = upper end
case <- c(rep("prior", length(cal_param_names)), rep("posterior", length(cal_param_names)))
pe <- rep(0, length(cal_param_names) * 2)
lend <- rep(0, length(cal_param_names) * 2)
uend <- rep(0, length(cal_param_names) * 2)

ggplot.data <- data.frame(par = par, case = case, pe = pe, lend = lend, uend = uend)
CalibPar <- read.xlsx("Inputs/MasterTable.xlsx", sheet = "CalibPar")
ggplot.data$pe[ggplot.data$case == "prior"] <- CalibPar$pe
ggplot.data$lend[ggplot.data$case == "prior"] <- CalibPar$pe - CalibPar$lower
ggplot.data$uend[ggplot.data$case == "prior"] <- CalibPar$upper - CalibPar$pe

ggplot.data$pe[ggplot.data$case == "posterior"] <- apply(calib.post, 2, median)
ggplot.data$lend[ggplot.data$case == "posterior"] <- apply(calib.post, 2, median) - apply(calib.post, 2, quantile, probs = 0.025)
ggplot.data$uend[ggplot.data$case == "posterior"] <- apply(calib.post, 2, quantile, probs = 0.975) - apply(calib.post, 2, median)

ggplot.data$case <- factor(ggplot.data$case, levels = c("prior", "posterior"))

p <- ggplot(ggplot.data, aes(x = case, y = pe, color = case)) +
  geom_point() +
  geom_errorbar(aes(ymin = pe - lend, ymax = pe + uend), width = .2, position = position_dodge(0.05)) +
  facet_wrap(~par, scales = "free_y", ncol = 5) +
  labs(y = "Value", x = "") +
  theme_bw()

## City-level validation, please go to CityLevelValidation.R ##
