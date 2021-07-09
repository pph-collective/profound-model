###############################################################################################
#######################        Model calibration runs         #################################
###############################################################################################
# Module for running model with uncalibrated data in parallel over calibration period
#
# Authors: Xiao Zang, PhD, Sam Bessey, MS
#
# People, Place and Health Collective, Department of Epidemiology, Brown University
#

rm(list=ls())

# ## Model setup parameters ##
seed         <- 2021    #define the initial seed, will draw different seeds based on the initial for calibration simulations
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  batch.ind <- 1
  outpath <- getwd()
  cores <- 1
} else{
  batch.ind <- strtoi(args[1])
  outpath <- strtoi(args[2])
  cores <- strtoi(args[3])
}
batch.size   <- 100000  #define the size of each batch of calibration simulations, default we have 10 batches, each with 100000 simulations

#Load required packages
print("Loading required packages")
library(dplyr)
library(openxlsx)
library(abind)
library(FME)
#Load packages for parallel performance
library(foreach)
library(doParallel)
# make and register the cluster given the provided number of cores
c1 = makeCluster(cores, outfile="")
registerDoParallel(c1)


#Load required scripts and functions
source("population.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("naloxone_available.R")
source("cost_effectiveness.R")
source("parallel.R")
source("prep_calibration_data.R")


# REVIEWED - check for both / is Calib_par_table.rds not used here? And do we know that if that file exists, the calibration sample data files necessarily exist?
## load or create calibration parameter sets for calibration simulation 
if(file.exists(paste0("Inputs/Calib_par_table.rds"))){
  #only load the indexed parameter set batch for calibration simulation
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds")) 
} else {
  ## Specify the number of calibration random parameter sets
  sample.size <- 1000000  #total number of calibration samples

  #load calibration parameter bounds and values
  CalibPar <- read.xlsx("Inputs/MasterTable.xlsx", "CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par

  # create and save calibration parameters using latin hypercube sampling with a set seed
  set.seed(5112021)
  calib.par <- Latinhyper(parRange, sample.size) %>% data.frame()
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))  #save sampled calibration parameter values

  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds"))
  rm(calib.par)
}

# generate stepwise seeds for calibration starting at initial seed
calib.seed.vt  <- seed + c(((batch.ind-1)*batch.size + 1):(batch.ind*batch.size))
# initialize calibratiobration_results table
calib.rs.table <- matrix(0, nrow = length(Calibration.data.ls), ncol = 15) 
colnames(calib.rs.table) <- c("index", "seed",
                              "od.death16", "od.death17", "od.death18", "od.death19",
                              "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                              "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19", 
                              "gof")
calib.rs.table[ , "index"] <- c(((batch.ind-1)*batch.size + 1):(batch.ind*batch.size))
calib.rs.table[ , "seed"]  <- calib.seed.vt

#initializbration_results matrix to savbration_results from parallel simulation
calibration_results <- matrix(0, nrow = length(Calibration.data.ls), ncol = 12)

v.rgn <- Calibration.data.ls[[1]]$v.rgn  #load vector for regions (required by simulation)

# parallel calibration simulation
# TODO: look into apply functions for parallel
calibration_results <- foreach(ss = 1:length(Calibration.data.ls), .combine = rbind, .packages= c('dplyr', 'abind')) %dopar% {
  yr_start    <- 2016   #simulation first year
  yr_end     <- 2020   #simulation last year
  ppl_info    <- c("sex", "race", "age", "residence", "curr.state",
                   "OU.state", "init.age", "init.state", "ever.od", "fx")            # information for each model individual
  agent_states     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
  v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
  num_states     <- length(agent_states)                                                     # number of states
  num_years        <- yr_end-yr_start+1
  timesteps         <- 12 * num_years                                                           # number of time cycles (in month)
  n.rgn       <- length(v.rgn)                                                       # number of regions
  
  # OUTPUT matrices and vectors
  v.od        <- rep(0, times = timesteps)                                                 # count of overdose events at each time step
  v.oddeath   <- rep(0, times = timesteps)                                                 # count of overdose deaths at each time step
  m.oddeath   <- matrix(0, nrow = timesteps, ncol = n.rgn)
  colnames(m.oddeath) <- v.rgn
  v.odpriv    <- rep(0, times = timesteps)                                                 # count of overdose events occurred at private setting at each time step
  v.odpubl    <- rep(0, times = timesteps)                                                 # count of overdose events occurred at public setting at each time step
  v.deathpriv <- rep(0, times = timesteps)                                                 # count of overdose deaths occurred at private setting at each time step
  v.deathpubl <- rep(0, times = timesteps)                                                 # count of overdose deaths occurred at public setting at each time step
  v.str       <- c("SQ", "Expand100")                                                # store the strategy names
  d.c         <- 0.03                                                                # discounting of costs by 3%
  cost.item   <- c("TotalCost", "NxCost")
  cost.matrix <- matrix(0, nrow=timesteps, ncol = length(cost.item))
  colnames(cost.matrix) <- cost.item
  m.oddeath.fx <- rep(0, times = timesteps)                                                # count of overdose deaths with fentanyl present at each time step
  m.oddeath.op <- rep(0, times = timesteps)                                                # count of overdose deaths among opioid users at each time step
  m.oddeath.st <- rep(0, times = timesteps)                                                # count of overdose deaths among stimulant users at each time step
  m.EDvisits   <- rep(0, times = timesteps)                                                # count of opioid overdose-related ED visits at each time step
  m.oddeath.hr <- rep(0, times = timesteps)                                                # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
  
  ## Initialize the study population - people who are at risk of opioid overdose
  ppl_info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
  init.pop  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
  
  outcomes <- parallel.fun(calib.seed = calib.seed.vt[ss], vparameters = Calibration.data.ls[[ss]])
}

calib.rs.table[,3:14] <- calibration_results   #pass calibratiobration_results tbration_results table
saveRDS(calib.rs.table, paste0("CalibrationOutputs", batch.ind, ".rds")) #save calibratiobration_results table to an rds, will combine all 10 tables/bacthes in a subsequent process

# stopCluster(c1)   #optional: stop clustering (breaking programs into different cores)