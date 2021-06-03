###############################################################################################
#######################           Model calibration           #################################
###############################################################################################
rm(list=ls())

# ## Model setup parameters ##
seed         <- 2021
batch.ind    <- 1       #specify the index of batch for calibration analysis
batch.size   <- 100000

#Load required packages
library(dplyr)
library(openxlsx)
library(abind)
#Load packages for parallel performance
library(foreach)
library(doParallel)
#Specify number of cores
# ncores = 5; c1 <- makeCluster(ncores)   #make clusters for local machine
c1 = Sys.getenv("SLURM_NTASKS")           #make clusters for clusters
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
if(file.exists(paste0("Inputs/Calib_par_table.rds"))){
  # calib.par           <- readRDS(paste0("Inputs/Calib_par_table.rds"))
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds"))
} else if (!file.exists(paste0("Inputs/Calib_par_table.rds"))){
  ## Specify the number of calibration random parameter sets
  sample.size <- 1000000
  batch.size  <- 100000
  library(FME)
  WB       <- loadWorkbook("Inputs/MasterTable.xlsx")
  CalibPar <- read.xlsx(WB, sheet="CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par
  set.seed(5112021)
  calib.par <- Latinhyper(parRange, sample.size)
  calib.par <- data.frame(calib.par)
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))
  source("Profound-CalibrationDataPrep.R")
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds"))
}

# rm(calib.par)

# nm.calp <- names(calib.par)
calib.seed.vt <- seed + c(((batch.ind-1)*batch.size + 1):(batch.ind*batch.size))

calib.rs.table <- matrix(0, nrow = length(Calibration.data.ls), ncol = 15)
colnames(calib.rs.table) <- c("index", "seed",
                              "od.death16", "od.death17", "od.death18", "od.death19",
                              "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                              "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19", 
                              "gof")
calib.rs.table[ , "index"] <- c(((batch.ind-1)*batch.size + 1):(batch.ind*batch.size))
calib.rs.table[ , "seed"] = calib.seed.vt

calib.results <- matrix(0, nrow = length(Calibration.data.ls), ncol = 12)
v.rgn <- Calibration.data.ls[[1]]$v.rgn
# export.par <- c("nlx.avail")

# calib.results <- foreach(ss = 1:nrow(calib.par), .combine = rbind, .errorhandling = 'remove', .packages= c('dplyr', 'abind')) %dopar% {
# calib.results <- foreach(ss = 1:10, .combine = rbind, .packages= c('dplyr', 'abind', 'openxlsx'), .export = ls(globalenv())) %dopar% {
# calib.results <- foreach(ss = 1:10, .combine = rbind, .packages= c('dplyr', 'abind', 'openxlsx'), .export = export.par) %dopar% {
calib.results <- foreach(ss = 1:length(Calibration.data.ls), .combine = rbind, .packages= c('dplyr', 'abind')) %dopar% {
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

calib.rs.table[3:14] <- calib.results
saveRDS(calib.rs.table, paste0("CalibrationOutputs", batch.ind, ".rds"))

stopCluster(c1)