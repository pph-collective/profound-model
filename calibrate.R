#!/usr/bin/env Rscript

###############################################################################################
#######################        Model calibration runs         #################################
###############################################################################################
###############################################################################################
####    1st step of Calibration: create random parameter samples and simulate over them    ####
####    3 sets of targets:  annual # overdose deaths:       od.deathYR (YR for year)       ####
####                        % of overdose death with fentanyl present:    fx.deathYR       ####
####                        annual # ED visits due to overdose:           ed.visitYR       ####
####    Calibrate multiple parameters (calib.par) simultaneously                           ####
####    Targets at the state level, 2016-2019                                              ####
####    Method: random calibration with Latin hypercube sampling                           ####
####    DUE to large samples, run calibration simulation with parallel in batches          ####
####    AFTER finishing all batches, run ResultAnalysis.R for the next step                ####
###############################################################################################
rm(list = ls())
# Module for running model with uncalibrated data in parallel over calibration period
#
# Authors: Xiao Zang, PhD, Sam Bessey, MS
#
# People, Place and Health Collective, Department of Epidemiology, Brown University
#

rm(list = ls())
print("Loading required packages")
library(dplyr)
library(openxlsx)
library(abind)
library(FME)
# Load packages for parallel performance
library(foreach)
library(doParallel)
library("argparser")

# ## Model setup parameters ##
args <- arg_parser("arguments")
args <- add_argument(args, "--seed", help = "seed for latin hypercube sampling", default = 5112021)
args <- add_argument(args, "--outfile", help = "file to store outputs", default = "OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--params", help = "name of the file with calibration params", default = "Inputs/CalibrationSampleData1.rds")
args <- add_argument(args, "--cores", help = "number of cores for parallel processing", default = 4)
args <- add_argument(args, "--index", help = "index of batch", default = 1)
args <- add_argument(args, "--size", help = "batch size", default = 1000)

argv <- parse_args(args)
seed <- as.integer(argv$seed)
cores <- as.integer(argv$cores)
batch.ind <- 4
batch.size <- as.integer(argv$size)
# note: outfile, params not used

# make and register the cluster given the provided number of cores
c1 <- makeCluster(cores, outfile = "")
registerDoParallel(c1)


# Load required scripts and functions
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")
source("parallel.R")

# REVIEWED - check for both / is Calib_par_table.rds not used here? And do we know that if that file exists, the calibration sample data files necessarily exist?
## load or create calibration parameter sets for calibration simulation
if (file.exists(paste0("calibration/prep_calibration_data/Calib_par_table.rds"))) {
  # only load the indexed parameter set batch for calibration simulation
  Calibration.data.ls <- readRDS(paste0("calibration/prep_calibration_data/CalibrationSampleData", batch.ind, ".rds"))
} else {
  ## Specify the number of calibration random parameter sets
  sample.size <- 10000 # total number of calibration samples
  source("prep_calibration_data.R")
  Calibration.data.ls <- readRDS(paste0("calibration/prep_calibration_data/CalibrationSampleData", batch.ind, ".rds"))
  # # load calibration parameter bounds and values
  # CalibPar <- read.xlsx("Inputs/MasterTable.xlsx", "CalibPar")
  # parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  # row.names(parRange) <- CalibPar$par
  # 
  # # create and save calibration parameters using latin hypercube sampling with a set seed
  # set.seed(seed)
  # calib.par <- Latinhyper(parRange, sample.size) %>% data.frame()
  # saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds")) # save sampled calibration parameter values
  # 
  # Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds"))
  
}

# generate stepwise seeds for calibration starting at initial seed
calib.seed.vt <- seed + c(((batch.ind - 1) * batch.size + 1):(batch.ind * batch.size))
# initialize calibratiobration_results table
calibration_results <- matrix(0, nrow = length(Calibration.data.ls), ncol = 15)
colnames(calibration_results) <- c(
  "index", "seed",
  "od.death16", "od.death17", "od.death18", "od.death19",
  "fx.death16", "fx.death17", "fx.death18", "fx.death19",
  "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19",
  "gof"
)
calibration_results[, "index"] <- c(((batch.ind - 1) * batch.size + 1):(batch.ind * batch.size))
calibration_results[, "seed"] <- calib.seed.vt

# initializbration_results matrix to savbration_results from parallel simulation
calibration_temp <- matrix(0, nrow = length(Calibration.data.ls), ncol = 12)
colnames(calibration_temp) <- c("od.death16", "od.death17", "od.death18", "od.death19",
                                   "fx.death16", "fx.death17", "fx.death18", "fx.death19",
                                   "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19")  
v.region <- Calibration.data.ls[[1]]$v.region # load vector for regions (required by simulation)

for (jj in 1:length(Calibration.data.ls)){
  Calibration.data.ls[[jj]]$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
}
# parallel calibration simulation
# TODO: look into apply functions for parallel

for (ss in 1:length(Calibration.data.ls)){
  print(paste0("Batch: ", batch.ind, "; simulation: ", ss))
  yr_start <- 2016 # simulation first year
  yr_end <- 2020 # simulation last year
  ppl_info <- c(
    "sex", "race", "age", "residence", "curr.state",
    "OU.state", "init.age", "init.state", "ever.od", "fx"
  ) # information for each model individual
  agent_states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead") # vector for state names
  v.oustate <- c("preb", "il.lr", "il.hr") # vector for active opioid use state names
  num_states <- length(agent_states) # number of states
  num_years <- yr_end - yr_start + 1
  timesteps <- 12 * num_years # number of time cycles (in month)
  num_regions <- length(v.region) # number of regions
  
  # OUTPUT matrices and vectors
  v.od <- rep(0, times = timesteps) # count of overdose events at each time step
  v.oddeath <- rep(0, times = timesteps) # count of overdose deaths at each time step
  m.oddeath <- matrix(0, nrow = timesteps, ncol = num_regions)
  colnames(m.oddeath) <- v.region
  v.odpriv <- rep(0, times = timesteps) # count of overdose events occurred at private setting at each time step
  v.odpubl <- rep(0, times = timesteps) # count of overdose events occurred at public setting at each time step
  v.deathpriv <- rep(0, times = timesteps) # count of overdose deaths occurred at private setting at each time step
  v.deathpubl <- rep(0, times = timesteps) # count of overdose deaths occurred at public setting at each time step
  v.nlxused <- rep(0, times = timesteps) # count of overdose deaths occurred at public setting at each time step
  v.str <- c("SQ", "Expand100") # store the strategy names
  d.c <- 0.03 # discounting of costs by 3%
  cost.item <- c("TotalCost", "NxCost")
  cost.matrix <- matrix(0, nrow = timesteps, ncol = length(cost.item))
  colnames(cost.matrix) <- cost.item
  m.oddeath.fx <- rep(0, times = timesteps) # count of overdose deaths with fentanyl present at each time step
  m.oddeath.op <- rep(0, times = timesteps) # count of overdose deaths among opioid users at each time step
  m.oddeath.st <- rep(0, times = timesteps) # count of overdose deaths among stimulant users at each time step
  m.EDvisits <- rep(0, times = timesteps) # count of opioid overdose-related ED visits at each time step
  m.oddeath.hr <- rep(0, times = timesteps) # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
  
  ## Initialize the study population - people who are at risk of opioid overdose
  ppl_info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
  init_ppl <- readRDS(paste0("Inputs/init_pop.rds"))
  
  outcomes <- parallel.fun(calib.seed = calib.seed.vt[ss], params = Calibration.data.ls[[ss]])
  calibration_temp[ss, ] <- outcomes
}
# calibration_temp <- foreach(ss = 1:length(Calibration.data.ls), .combine = rbind, .packages = c("dplyr", "abind"), .errorhandling = "remove") %dopar% {
#   yr_start <- 2016 # simulation first year
#   yr_end <- 2020 # simulation last year
#   ppl_info <- c(
#     "sex", "race", "age", "residence", "curr.state",
#     "OU.state", "init.age", "init.state", "ever.od", "fx"
#   ) # information for each model individual
#   agent_states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead") # vector for state names
#   v.oustate <- c("preb", "il.lr", "il.hr") # vector for active opioid use state names
#   num_states <- length(agent_states) # number of states
#   num_years <- yr_end - yr_start + 1
#   timesteps <- 12 * num_years # number of time cycles (in month)
#   num_regions <- length(v.region) # number of regions
# 
#   # OUTPUT matrices and vectors
#   v.od <- rep(0, times = timesteps) # count of overdose events at each time step
#   v.oddeath <- rep(0, times = timesteps) # count of overdose deaths at each time step
#   m.oddeath <- matrix(0, nrow = timesteps, ncol = num_regions)
#   colnames(m.oddeath) <- v.region
#   v.odpriv <- rep(0, times = timesteps) # count of overdose events occurred at private setting at each time step
#   v.odpubl <- rep(0, times = timesteps) # count of overdose events occurred at public setting at each time step
#   v.deathpriv <- rep(0, times = timesteps) # count of overdose deaths occurred at private setting at each time step
#   v.deathpubl <- rep(0, times = timesteps) # count of overdose deaths occurred at public setting at each time step
#   v.nlxused <- rep(0, times = timesteps) # count of overdose deaths occurred at public setting at each time step
#   v.str <- c("SQ", "Expand100") # store the strategy names
#   d.c <- 0.03 # discounting of costs by 3%
#   cost.item <- c("TotalCost", "NxCost")
#   cost.matrix <- matrix(0, nrow = timesteps, ncol = length(cost.item))
#   colnames(cost.matrix) <- cost.item
#   m.oddeath.fx <- rep(0, times = timesteps) # count of overdose deaths with fentanyl present at each time step
#   m.oddeath.op <- rep(0, times = timesteps) # count of overdose deaths among opioid users at each time step
#   m.oddeath.st <- rep(0, times = timesteps) # count of overdose deaths among stimulant users at each time step
#   m.EDvisits <- rep(0, times = timesteps) # count of opioid overdose-related ED visits at each time step
#   m.oddeath.hr <- rep(0, times = timesteps) # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
# 
#   ## Initialize the study population - people who are at risk of opioid overdose
#   ppl_info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
#   init_ppl <- readRDS(paste0("Inputs/init_pop.rds"))
# 
#   outcomes <- parallel.fun(calib.seed = calib.seed.vt[ss], params = Calibration.data.ls[[ss]])
# }

calibration_results[, 3:14] <- calibration_temp #  calibration results
saveRDS(calibration_results, paste0("calibration/CalibrationOutputs", batch.ind, ".rds")) # save calibration_results table to an rds, will combine all 10 tables/bacthes in a subsequent process

# stopCluster(c1)   #optional: stop clustering
