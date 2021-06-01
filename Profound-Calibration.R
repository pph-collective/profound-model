###############################################################################################
#######################           Model calibration           #################################
###############################################################################################
rm(list=ls())

# ## Model setup parameters ##
seed         <- 2021

#Load required packages
library(dplyr)
library(openxlsx)
library(abind)
#Load packages for parallel performance
library(foreach)
library(doParallel)
#Specify number of cores
ncores = 5
c1 <- makeCluster(ncores)
# c1 = Sys.getenv("SLURM_NTASKS")
registerDoParallel(c1)


#Load scripts
source("Profound-Function-PopInitialization.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-Function-Microsimulation.R")
source("Profound-DecisionTree.R")
source("Profound-Function-NxAvailAlgm.R")
source("Profound-CEA.R")
source("Profound-Function-Parallel.R")

## Model parameter updates for calibration process ##
if(file.exists(paste0("Inputs/Calib_par_table.rds")) & file.exists(paste0("Inputs/CalibrationSampleData.rds"))){
  calib.par           <- readRDS(paste0("Inputs/Calib_par_table.rds"))
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData.rds"))
} else if (!file.exists(paste0("Inputs/Calib_par_table.rds"))){
  library(FME)
  CalibPar <- read.xlsx(WB, sheet="CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par
  set.seed(5112021)
  calib.par <- Latinhyper(parRange, 10000)
  calib.par <- data.frame(calib.par)
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))
  source("Profound-CalibrationDataPrep.R")
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData.rds"))
}

nm.calp <- names(calib.par)
calib.seed.vt <- seed + c(1:nrow(calib.par))

calib.rs.table <- matrix(0, nrow = nrow(calib.par), ncol = 15+length(nm.calp))
colnames(calib.rs.table) <- c("index", "seed", nm.calp, 
                              "od.death16", "od.death17", "od.death18", "od.death19",
                              "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                              "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19", "gof")
calib.rs.table[ , "index"] <- c(1:nrow(calib.par))
calib.rs.table[ , 3:(2+length(nm.calp))] <- as.matrix(calib.par)
calib.rs.table[ , "seed"] = calib.seed.vt

calib.results <- matrix(0, nrow = 10, ncol = 12)
v.rgn <- Calibration.data.ls[[1]]$v.rgn
export.par <- c("nlx.avail")

# calib.results <- foreach(ss = 1:nrow(calib.par), .combine = rbind, .errorhandling = 'remove', .packages= c('dplyr', 'abind')) %dopar% {
# calib.results <- foreach(ss = 1:10, .combine = rbind, .packages= c('dplyr', 'abind', 'openxlsx'), .export = ls(globalenv())) %dopar% {
calib.results <- foreach(ss = 1:10, .combine = rbind, .packages= c('dplyr', 'abind', 'openxlsx'), .export = export.par) %dopar% {
# calib.results <- foreach(ss = 1:10, .combine = rbind, .packages= c('dplyr', 'abind')) %dopar% {
  yr.first    <- 2016
  yr.last     <- 2020
  pop.info    <- c("sex", "race", "age", "residence", "curr.state",
                   "OU.state", "init.age", "init.state", "ever.od", "fx")            # information for each model individual
  v.state     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
  v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
  n.state     <- length(v.state)                                                     # number of states
  n.yr        <- yr.last-yr.first+1
  n.t         <- 12 * n.yr                                                           # number of time cycles (in month)
  n.rgn       <- length(v.rgn)                                                       # number of regions
  
  # OUTPUT matrices and vectors
  v.od        <- rep(0, times = n.t)                                                 # count of overdose events at each time step
  v.oddeath   <- rep(0, times = n.t)                                                 # count of overdose deaths at each time step
  m.oddeath   <- matrix(0, nrow = n.t, ncol = n.rgn)
  colnames(m.oddeath) <- v.rgn
  v.odpriv    <- rep(0, times = n.t)                                                 # count of overdose events occurred at private setting at each time step
  v.odpubl    <- rep(0, times = n.t)                                                 # count of overdose events occurred at public setting at each time step
  v.deathpriv <- rep(0, times = n.t)                                                 # count of overdose deaths occurred at private setting at each time step
  v.deathpubl <- rep(0, times = n.t)                                                 # count of overdose deaths occurred at public setting at each time step
  v.str       <- c("SQ", "Expand100")                                                # store the strategy names
  d.c         <- 0.03                                                                # discounting of costs by 3%
  cost.item   <- c("TotalCost", "NxCost")
  cost.matrix <- matrix(0, nrow=n.t, ncol = length(cost.item))
  colnames(cost.matrix) <- cost.item
  m.oddeath.fx <- rep(0, times = n.t)                                                # count of overdose deaths with fentanyl present at each time step
  m.oddeath.op <- rep(0, times = n.t)                                                # count of overdose deaths among opioid users at each time step
  m.oddeath.st <- rep(0, times = n.t)                                                # count of overdose deaths among stimulant users at each time step
  m.EDvisits   <- rep(0, times = n.t)                                                # count of opioid overdose-related ED visits at each time step
  m.oddeath.hr <- rep(0, times = n.t)                                                # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
  
  ## Initialize the study population - people who are at risk of opioid overdose
  pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
  init.pop  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
  
  outcomes <- parallel.fun(calib.seed = calib.seed.vt[ss], vparameters = Calibration.data.ls[[ss]])
  outcomes
}

stopCluster(c1)




## Use a regular FOR loop ##
for (ss in 1:nrow(calib.par)){
  print(paste0("Parameter set: ", ss))
  calib.seed = calib.rs.table[ss , "seed"]
  for (pp in 1:length(nm.calp)){
    assign(nm.calp[pp], calib.par[ss,pp])
  }
  
  # Update Overdose probability matrix (per month)
  od.matrix             <- matrix(0, nrow = 4, ncol = 2)
  rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
  colnames(od.matrix)   <- c("first", "subs")
  od.matrix["preb", "subs"]   <- od.preb.sub
  od.matrix["il.lr", "subs"]  <- od.il.lr.sub
  od.matrix["il.hr", "subs"]  <- od.il.lr.sub * multi.hr
  od.matrix["NODU", "subs"]   <- od.NODU.sub
  od.matrix[ , "first"]       <- od.matrix[ , "subs"] / multi.sub
  # Update decisiom tree parameters
  OD_911_pub    <- OD_911_priv * OD_911_pub_mul
  
  ## Run the simulation model over all sets of parameters
  sim_sq    <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = calib.seed)        # run for status quo
  
  calib.rs.table[ss, c("od.death16", "od.death17", "od.death18", "od.death19")] <- colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5))[-5]
  calib.rs.table[ss, c("fx.death16", "fx.death17", "fx.death18", "fx.death19")] <- (colSums(matrix(sim_sq$m.oddeath.fx, nrow=12, ncol=5)) / colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5)))[-5]
  calib.rs.table[ss, c("ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")] <- colSums(matrix(sim_sq$m.EDvisits, nrow=12, ncol=5))[-5]
}

Target   <- read.xlsx(WB, sheet="Target")
tar.data <- Target$pe

for (ss in 1:nrow(calib.rs.table)){
  model.pred <- calib.rs.table[ss, c("od.death16", "od.death17", "od.death18", "od.death19",
                                     "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                                     "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")]
  gof <- 0
  for (j in 1:length(tar.data)){
    if (j %in% c(5:8)){
      gof <- gof + (abs(model.pred[j]*100 - tar.data[j]*100)/(tar.data[j]*100))/length(tar.data)
    } else {
      gof <- gof + (abs(model.pred[j] - tar.data[j])/tar.data[j])/length(tar.data)
    }
  }
  calib.rs.table[ss, "gof"] <- gof
}

sorted.mx <- calib.rs.table[order(calib.rs.table[ , "gof"], decreasing = F), ]

library(xlsx)
write.xlsx(sorted.mx, file="Calibration_preliminary.xlsx", 
           col.names = T, row.names = F)

calib.sp <- 20
calib.result.mx <- sorted.mx[1:calib.sp, ]

par(mfrow=c(1,3))

md.oddeath <- apply(calib.result.mx[ , c("od.death16", "od.death17", "od.death18", "od.death19")], 2, median)
ymax       <- max(calib.result.mx[ , c("od.death16", "od.death17", "od.death18", "od.death19")])
plot(x = 2016:2019,  tar.data[1:4], col = "black" , pch = 18, xlab="Year", ylab="Number of opioid overdose deaths", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', ylim = c(0, 500), frame.plot = FALSE)
# title("A", adj = 0, line =-0.5)
for (i in 1:calib.sp){
  lines(x = 2016:2019, y = calib.result.mx[i, c("od.death16", "od.death17", "od.death18", "od.death19")], col = adjustcolor("indianred1", alpha = 0.2), lwd=2)
}
lines(x = 2016:2019, y = md.oddeath, col = "red", lwd=3)
points(x = 2016:2019,  tar.data[1:4], col = "black" , pch = 16, cex=1.2, cex.axis=0.95)
axis(1, at = 2016:2019, pos=0, lwd.ticks=0, cex.axis=1.2)
# axis(side=1, pos=0, lwd.ticks=0)
abline(h=0)
# legend(x="bottomright",
#        col=c("dodgerblue","firebrick2", "lightgrey"), 
#        lwd=c(1.2, 1, 0.5), lty = c(1, NA, NA), pch=c(16, 16,16), pt.cex = c(1,1.1,1), cex = 0.8, bty = "n")


md.fxoddeath <- apply(calib.result.mx[ ,c("fx.death16", "fx.death17", "fx.death18", "fx.death19")], 2, median)
ymax       <- max(calib.result.mx[ , c("fx.death16", "fx.death17", "fx.death18", "fx.death19")])
plot(x = 2016:2019,  tar.data[5:8], col = "black" , pch = 18, xlab="Year", ylab="Proportion of OOD deaths with fentanyl present", cex=1.2, cex.axis=1.2, cex.lab = 1.3,
     xaxt='n', yaxt='n', ylim = c(0, 1), frame.plot = FALSE)
# title("A", adj = 0, line =-0.5)
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
# title("A", adj = 0, line =-0.5)
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
legend("top", 
       legend = c("Target", "Model: median", "Model: simulation"), 
       col = c("black",  "red", adjustcolor("indianred1", alpha = 0.2)), 
       pch = c(16, NA, NA), 
       lty = c(NA, 1, 1), 
       lwd=3,
       bty = "n", 
       pt.cex = 2, 
       cex = 1.2, 
       text.col = "black", 
       horiz = T)


library(ggplot2)
calib.post  <- calib.result.mx[ , 3:(2+length(nm.calp))]
par  <- rep(colnames(calib.post), 2)
case <- c(rep("prior", length(nm.calp)), rep("posterior", length(nm.calp)))
pe   <- rep(0, length(nm.calp)*2)
lend <- rep(0, length(nm.calp)*2)
uend <- rep(0, length(nm.calp)*2)
ggplot.data <- data.frame(par = par, case = case, pe = pe, lend = lend, uend = uend)
CalibPar <- read.xlsx(file = "Inputs/MasterTable.xlsx", sheetIndex = 20)
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
