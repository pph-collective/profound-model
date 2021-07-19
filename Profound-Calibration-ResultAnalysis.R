###########################################################################################################
#######################                Model calibration - Result analysis         ########################
###########################################################################################################

###########################################################################################################
####    2nd step of calibration: Analyze calibration simulation results and select best-fit subsets    ####
####    3 sets of targets:  annual # overdose deaths:       od.deathYR (YR for year)                   ####
####                        % of overdose death with fentanyl present:    fx.deathYR                   ####
####                        annual # ED visits due to overdose:           ed.visitYR                   ####
####    Calibrate multiple parameters (calib.par) simultaneously                                       ####
####    Targets at the state level, 2016-2019                                                          ####
####    Additional validation for OD deaths at the city level at the bottom                            ####
####    Selection criteria: goodness-of-fit (mean percentage deviation)                                ####
####    One target, % fentanyl presence in od deaths, is a %, will *100 to ensure comparability        ####
####    Function to plot calibration results against targets also included                             ####
###########################################################################################################

rm(list=ls())
#Load required packages
library(dplyr)
library(openxlsx)
library(abind)

#Load scripts
source("Profound-Function-PopInitialization.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-Function-Microsimulation.R")
source("Profound-DecisionTree.R")
source("Profound-Function-NxAvailAlgm.R")
source("Profound-CEA.R")
source("Profound-Function-Parallel.R")

calib.rs.table <- NULL             #initiate calibration result table as an empty table

# Load calibration results from 10 batches of simulations (from previous step)
for (batch.ind in 1:10){
  calib.rs.temp  <- readRDS(paste0("calibration/CalibrationOutputs", batch.ind, ".rds"))
  calib.rs.table <- rbind(calib.rs.table, calib.rs.temp)
}
rm(calib.rs.temp)

# Load table for calibration parameters (randomly drawn values from Latin Hypercube sampling from last step)
calib.par <- readRDS(paste0("calibration/Calib_par_table.rds"))
# Combine calibration results and parameter values
calib.rs.table <- cbind(calib.rs.table, calib.par)
# Load calibration target data
WB       <- loadWorkbook("Inputs/MasterTable.xlsx")
Target   <- read.xlsx(WB, sheet="Target")
tar.data <- Target$pe

# Calculate goodness of fit (gof) for each calibration simulation, gof based on mean percenatge deviation
for (ss in 1:nrow(calib.rs.table)){
  model.pred <- calib.rs.table[ss, c("od.death16", "od.death17", "od.death18", "od.death19",
                                     "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                                     "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")]
  gof <- 0
  for (j in 1:length(tar.data)){
    if (j %in% c(5:8)){   # targets 5-8 are % measures, % of overdose deaths with fentanyl present, times that with 100 to make it comparable with other targets
      gof <- gof + (abs(model.pred[j]*100 - tar.data[j]*100)/(tar.data[j]*100))/length(tar.data)
    } else {
      gof <- gof + (abs(model.pred[j] - tar.data[j])/tar.data[j])/length(tar.data)
    }
  }
  calib.rs.table[ss, "gof"] <- gof
}

# Sort calibration results table according to god (ascending order)
sorted.mx <- calib.rs.table[order(calib.rs.table[ , "gof"], decreasing = F), ]

# Subselect target number of best-fit subsets and save them in excel
library(xlsx)
calib.sp <- 500    #define target number of best-fitting subsets from calibration for subsequent analysis
calib.result.mx <- sorted.mx[1:calib.sp, ]
write.xlsx(calib.result.mx, file=paste0("calibration/Calibrated_results.xlsx"),
           col.names = T, row.names = F)


## Save calibrated results as parameter lists in 'vparameters' (combine with uncalibrated parameters) ##
#INPUT PARAMETERS
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"
source("Profound-DataInput.R")
nm.calp <- names(calib.par)
calibrated.parameters <- list()

for (cc in 1:nrow(calib.result.mx)){
  for (pp in 1:length(nm.calp)){
    vparameters[[nm.calp[pp]]] <- calib.result.mx[cc, nm.calp[pp]]
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
  
  calibrated.parameters[[cc]] <- vparameters
}
calibrated.seed <- calib.result.mx[,"seed"]
saveRDS(calibrated.parameters, file = paste0("calibration/CalibratedData.rds"))
saveRDS(calibrated.seed, file = paste0("calibration/CalibratedSeed.rds"))


##############################################
## Plot calibration targets against targets ##
##############################################
par(mfrow=c(1,3))

md.oddeath <- apply(calib.result.mx[ , c("od.death16", "od.death17", "od.death18", "od.death19")], 2, mean)
ymax       <- max(calib.result.mx[ , c("od.death16", "od.death17", "od.death18", "od.death19")])
plot(x = 2016:2019,  tar.data[1:4], col = "black" , pch = 18, xlab="Year", ylab="Number of opioid overdose deaths", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', ylim = c(0, 500), frame.plot = FALSE)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("od.death16", "od.death17", "od.death18", "od.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=3)
}
lines(x = 2016:2019, y = md.oddeath, col = "red", lwd=4)
points(x = 2016:2019,  tar.data[1:4], col = "black" , pch = 16, cex=1.2, cex.axis=0.95)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
abline(h=0)
legend("top",
       legend = c("Target", "Model: mean", "Model: simulation"),
       col = c("black",  "red", adjustcolor("indianred1", alpha = 0.2)),
       pch = c(16, NA, NA),
       lty = c(NA, 1, 1),
       lwd=3,
       bty = "n",
       pt.cex = 1.2,
       cex = 1.1,
       text.col = "black",
       horiz = T)


md.fxoddeath <- apply(calib.result.mx[ ,c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], 2, median)
ymax       <- max(calib.result.mx[ , c("fx.death16", "fx.death17", "fx.death18", "fx.death19")])
plot(x = 2016:2019,  tar.data[5:8], col = "black" , pch = 18, xlab="Year", ylab="Proportion of OOD deaths with fentanyl present", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', yaxt='n', ylim = c(0, 1), frame.plot = FALSE)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=2)
}
lines(x = 2016:2019, y = md.fxoddeath, col = "red", lwd=3)
points(x = 2016:2019,  tar.data[5:8], col = "black" , pch = 16, cex=1.2, cex.axis=0.95)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c("0", "20%", "40%", "60%", "80%", "100%"), cex.axis=1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)


md.edvisits <- apply(calib.result.mx[ ,c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], 2, median)
ymax       <- max(calib.result.mx[ , c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")])
plot(x = 2016:2019,  tar.data[9:12], col = "black" , pch = 18, xlab="Year", ylab="Number of ED visists for opioid overdose", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', ylim = c(0, 2500), frame.plot = FALSE)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=2)
}
lines(x = 2016:2019, y = md.edvisits, col = "red", lwd=3)
points(x = 2016:2019,  tar.data[9:12], col = "black" , pch = 16, cex=1.2, cex.axis=1.2)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
# axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c("0", "20%", "40%", "60%", "80%", "100%"))
# axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)

par(mfrow=c(1,1))
legend("bottom", 
       legend = c("Target", "Model: mean", "Model: simulation"), 
       col = c("black",  "red", adjustcolor("indianred1", alpha = 0.2)), 
       pch = c(16, NA, NA), 
       lty = c(NA, 1, 1), 
       lwd=3,
       bty = "n", 
       pt.cex = 2, 
       cex = 1.1, 
       text.col = "black", 
       horiz = T)


##############################################################################
#### Plot model parameter values prior and posterior to model calibration ####
##############################################################################
nm.calp <- names(calib.par)
library(ggplot2)
calib.post  <- calib.result.mx[ , (dim(calib.result.mx)[2]-length(nm.calp)+1):dim(calib.result.mx)[2]]
par  <- rep(colnames(calib.post), 2)
case <- c(rep("prior", length(nm.calp)), rep("posterior", length(nm.calp)))
pe   <- rep(0, length(nm.calp)*2)
lend <- rep(0, length(nm.calp)*2)
uend <- rep(0, length(nm.calp)*2)
ggplot.data <- data.frame(par = par, case = case, pe = pe, lend = lend, uend = uend)
CalibPar <- read.xlsx("Inputs/MasterTable.xlsx", sheet = "CalibPar")
ggplot.data$pe[ggplot.data$case == "prior"]   <- CalibPar$pe
ggplot.data$lend[ggplot.data$case == "prior"] <- CalibPar$pe - CalibPar$lower
ggplot.data$uend[ggplot.data$case == "prior"] <- CalibPar$upper - CalibPar$pe

ggplot.data$pe[ggplot.data$case == "posterior"]   <- apply(calib.post, 2, median)
ggplot.data$lend[ggplot.data$case == "posterior"] <- apply(calib.post, 2, median) - apply(calib.post, 2, quantile, probs = 0.025)
ggplot.data$uend[ggplot.data$case == "posterior"] <- apply(calib.post, 2, quantile, probs = 0.975) - apply(calib.post, 2, median)

ggplot.data$case <- factor(ggplot.data$case, levels = c("prior", "posterior"))

p<- ggplot(ggplot.data, aes(x=case, y=pe, color=case)) + 
  geom_point()+
  geom_errorbar(aes(ymin=pe-lend, ymax=pe+uend), width=.2, position=position_dodge(0.05)) +
  facet_wrap(~ par, scales = "free_y", ncol = 5) +
  labs(y="Value", x="") + theme_bw()

## FOR City-level validation, please go to CityLevelValidation.R ##