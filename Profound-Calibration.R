###############################################################################################
#######################           Model calibration           #################################
###############################################################################################
rm(list=ls())

# ## Model setup parameters ##
seed         <- 2021    #define the initial seed, will draw different seeds based on the initial for calibration simulations
# batch.ind    <- 1       #specify the index of batch/job for calibration analysis (TO SAM: this is where you may need to modify, from 1-10 when submitting each of the 10 batches/jobs)
args <- commandArgs(trailingOnly=TRUE)
if (length(args) == 0){
  cores <- 1
  batch.ind <- 1
  outpath <- getwd()
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
#Load packages for parallel performance
library(foreach)
library(doParallel)
#Specify number of cores
# ncores = 5; c1 <- makeCluster(ncores)   #make clusters for local machine
c1 = makeCluster(cores)           #make clusters for clusters (TO SAM: not sure if this is the right command, used for Compute Canada)
registerDoParallel(c1)                    #register cluster


#Load required scripts and functions
source("population.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("naloxone_available.R")
source("cost_effectiveness.R")
source("parallel.R")

## load or create calibration parameter sets for calibration simulation 
#(TO SAM: all rds files were saved in Google Drive, may need to update the path, i.e. Inputs)
print("loading calibration parameters")
if(file.exists(paste0("Inputs/Calib_par_table.rds"))){
  #only load the indexed parameter set batch for calibration simulation
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds")) 
} else if (!file.exists(paste0("Inputs/Calib_par_table.rds"))){
  ## Specify the number of calibration random parameter sets
  sample.size <- 1000000  #total number of calibration samples
  batch.size  <- 100000   #number of samples per calibration batch
  library(FME)            #load package for latin hypercube function
  #load calibration parameter bounds and values
  WB       <- loadWorkbook("Inputs/MasterTable.xlsx")
  CalibPar <- read.xlsx(WB, sheet="CalibPar")
  parRange <- data.frame(min = CalibPar$lower, max = CalibPar$upper)
  row.names(parRange) <- CalibPar$par
  set.seed(5112021)
  calib.par <- Latinhyper(parRange, sample.size)            #use latin hypercube to draw random samples for parameters
  calib.par <- data.frame(calib.par)
  saveRDS(calib.par, paste0("Inputs/Calib_par_table.rds"))  #save sampled calibration parameter values
  source("prep_calibration_data.R")                  #prepare calibration data (as lists) and save them in rds files
  Calibration.data.ls <- readRDS(paste0("Inputs/CalibrationSampleData", batch.ind, ".rds"))
  rm(calib.par)
}

#generate seeds for calibration (incremental by 1 from the initial seed)
calib.seed.vt  <- seed + c(((batch.ind-1)*batch.size + 1):(batch.ind*batch.size))
#initialize calibration results table
calib.rs.table <- matrix(0, nrow = length(Calibration.data.ls), ncol = 15) 
colnames(calib.rs.table) <- c("index", "seed",
                              "od.death16", "od.death17", "od.death18", "od.death19",
                              "fx.death16", "fx.death17", "fx.death18", "fx.death19", 
                              "ed.visit16", "ed.visit17", "ed.visit18", "ed.visit19", 
                              "gof")
calib.rs.table[ , "index"] <- c(((batch.ind-1)*batch.size + 1):(batch.ind*batch.size))
calib.rs.table[ , "seed"]  <- calib.seed.vt

#initialize results matrix to save results from parallel simulation
calib.results <- matrix(0, nrow = length(Calibration.data.ls), ncol = 12)

v.rgn <- Calibration.data.ls[[1]]$v.rgn  #load vector for regions (required by simulation)

#parallel calibration simulation
calib.results <- foreach(ss = 1:length(Calibration.data.ls), .combine = rbind, .packages= c('dplyr', 'abind')) %dopar% {
  yr.first    <- 2016        # starting year of simulation 
  yr.last     <- 2020        # end year of simulation (also the year for evaluation)
  d.c         <- 0.03        # discounting of costs by 3%
  
  source("io_setup.R")
  
  ## Initialize the study population - people who are at risk of opioid overdose
  pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
  init_ppl  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
  
  outcomes <- parallel.fun(calib.seed = calib.seed.vt[ss], params = Calibration.data.ls[[ss]])
  outcomes
}

calib.rs.table[,3:14] <- calib.results   #pass calibration results to results table
saveRDS(calib.rs.table, paste0("CalibrationOutputs", batch.ind, ".rds")) #save calibration results table to an rds, will combine all 10 tables/bacthes in a subsequent process

# stopCluster(c1)   #optional: stop clustering (breaking programs into different cores)