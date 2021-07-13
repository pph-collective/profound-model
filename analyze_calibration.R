###############################################################################################
#######################         Calibration Analysis          #################################
###############################################################################################
# Module for analyzing calibration runs and selecting those with best goodness of fit
#
# Authors: Xiao Zang, PhD, Sam Bessey, MS
#
# People, Place and Health Collective, Department of Epidemiology, Brown University
#
###############################################################################################

rm(list = ls())
# Load required packages
library(dplyr)
library(openxlsx)
library(abind)
library(ggplot2)
library(xlsx)

# Load scripts
source("population.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("naloxone_available.R")
source("cost_effectiveness.R")
source("parallel.R")
source("data_input.R")


calibration_results <- NULL

for (batch.ind in 1:10) {
  temp_results <- readRDS(paste0("calibration/CalibrationOutputs", batch.ind, ".rds"))
  calibration_results <- rbind(calibration_results, temp_results)
}
rm(temp_results)

calibration_params <- readRDS(paste0("calibration/Calib_par_table.rds"))

calibration_results <- cbind(calibration_results, calibration_params)

# read in workbook
WB <- loadWorkbook("Inputs/MasterTable.xlsx")
Target <- read.xlsx(WB, sheet = "Target")
tar.data <- Target$pe

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

sorted.mx <- calibration_results[order(calibration_results[, "gof"], decreasing = F), ]

# REVIEWED sp=sample, choosing top samples. calib.sp. We seem to also use "sp" to refer to region-specific
calib.sp <- 100
calib.result.mx <- sorted.mx[1:calib.sp, ]
write.xlsx(calib.result.mx,
  file = paste0("calibration/Calibrated_results.xlsx"),
  col.names = T, row.names = F
)


## save calibrated results as parameter lists (prepare for main analysis)##
# INPUT PARAMETERS
# REVIEWED change ov and sp to "overall" and "regional"? sw.EMS.ODloc?
sw.EMS.ODloc <- "ov" # Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"
# REVIEWED nm.calp = parameter_names?
nm.calp <- names(calibration_params)
calibrated.parameters <- list()


# REVIEWED: random sampling, save list of parameters
# for each selected calibration run, determine the run's parameters
for (cc in 1:nrow(calib.result.mx)) {
  for (pp in 1:length(nm.calp)) {
    params[[nm.calp[pp]]] <- calib.result.mx[cc, nm.calp[pp]]
  }
  # Overdose probability parameter matrix (per month)
  od.matrix <- matrix(0, nrow = 4, ncol = 2)
  rownames(od.matrix) <- c("preb", "il.lr", "il.hr", "NODU")
  colnames(od.matrix) <- c("first", "subs")
  od.matrix["preb", "subs"] <- params$od.preb.sub
  od.matrix["il.lr", "subs"] <- params$od.il.lr.sub
  od.matrix["il.hr", "subs"] <- params$od.il.lr.sub * params$multi.hr
  od.matrix["NODU", "subs"] <- params$od.NODU.sub
  od.matrix[, "first"] <- od.matrix[, "subs"] / params$multi.sub
  params$od.matrix <- od.matrix

  # Baseline mortality parameters, excluding overdose (per month)
  mor.matrix <- matrix(0, nrow = 2, ncol = length(mor.gp))
  rownames(mor.matrix) <- c("bg", "drug")
  colnames(mor.matrix) <- mor.gp
  mor.matrix["bg", ] <- params$mor.bg
  mor.matrix["drug", ] <- params$mor.drug
  params$mor.matrix <- mor.matrix
  params$OD_911_pub <- params$OD_911_priv * params$OD_911_pub_mul

  calibrated.parameters[[cc]] <- params
}

# alibrated.seed <- calib.result.mx[,"seed"]
saveRDS(calibrated.parameters, file = paste0("calibration/CalibratedData.rds"))
saveRDS(calib.result.mx[, "seed"], file = paste0("calibration/CalibratedSeed.rds"))



## Plot ##
par(mfrow = c(1, 3))
# REVIEWED md = median, change to mean
mean_oddeath <- apply(calib.result.mx[, c("od.death16", "od.death17", "od.death18", "od.death19")], 2, mean)
ymax <- max(calib.result.mx[, c("od.death16", "od.death17", "od.death18", "od.death19")])
plot(
  x = 2016:2019, tar.data[1:4], col = "black", pch = 18, xlab = "Year", ylab = "Number of opioid overdose deaths", cex = 1.2, cex.axis = 1.2, cex.lab = 1.3,
  xaxt = "n", ylim = c(0, 500), frame.plot = FALSE
)
# title("A", adj = 0, line =-0.5)
for (i in 1:calib.sp) {
  lines(x = 2016:2019, y = calib.result.mx[i, c("od.death16", "od.death17", "od.death18", "od.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd = 3)
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


md.fxoddeath <- apply(calib.result.mx[, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], 2, median)
# REVIEWED: change from hardcoded ylim; mean instead of median for everything
plot(
  x = 2016:2019, tar.data[5:8], col = "black", pch = 18, xlab = "Year", ylab = "Proportion of OOD deaths with fentanyl present", cex = 1.2, cex.axis = 1.2, cex.lab = 1.3,
  xaxt = "n", yaxt = "n", ylim = c(0, 1), frame.plot = FALSE
)
# title("A", adj = 0, line =-0.5)
for (i in 1:calib.sp) {
  lines(x = 2016:2019, y = calib.result.mx[i, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd = 2)
}
lines(x = 2016:2019, y = md.fxoddeath, col = "red", lwd = 3)
points(x = 2016:2019, tar.data[5:8], col = "black", pch = 16, cex = 1.2, cex.axis = 0.95)
axis(1, at = 2016:2019, pos = 0, lwd.ticks = 0, cex.axis = 1.2)
axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8, 1), labels = c("0", "20%", "40%", "60%", "80%", "100%"), cex.axis = 1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h = 0)


md.edvisits <- apply(calib.result.mx[, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], 2, median)
plot(
  x = 2016:2019, tar.data[9:12], col = "black", pch = 18, xlab = "Year", ylab = "Number of ED visists for opioid overdose", cex = 1.2, cex.axis = 1.2, cex.lab = 1.3,
  xaxt = "n", ylim = c(0, 2500), frame.plot = FALSE
)

for (i in 1:calib.sp) {
  lines(x = 2016:2019, y = calib.result.mx[i, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], col = adjustcolor("indianred1", alpha = 0.2), lwd = 2)
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


nm.calp <- names(calibration_params)


calib.post <- calib.result.mx[, (dim(calib.result.mx)[2] - length(nm.calp) + 1):dim(calib.result.mx)[2]]
# REVIEWED par -> param
par <- rep(colnames(calib.post), 2)
# REVIEWED pe = point estimate, lend = lower end, uend = upper end
case <- c(rep("prior", length(nm.calp)), rep("posterior", length(nm.calp)))
pe <- rep(0, length(nm.calp) * 2)
lend <- rep(0, length(nm.calp) * 2)
uend <- rep(0, length(nm.calp) * 2)

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
