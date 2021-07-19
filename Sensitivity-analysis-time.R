#######################################################################################################
#####################            Sensitivity analysis - time horizon             ######################
#######################################################################################################
# To evaluate the sensitivity of program evaluation results in response to change in evaluation time horizon
# Simply by setting different "yr.last"

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
yr.first     <- 2016   #first year of evaluation
yr.last      <- 2024   #last year of evaluation
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

# Program data
library(openxlsx)
pg.data   <- read.xlsx("Ignore/OEND_program.xlsx", sheet = "Project Weber")
pg.levels <- c(1, 5, 10, 20, 50)
pg.add.array <- array(0, dim = c(length(pg.levels), 2, dim(pg.data)[1]))
dimnames(pg.add.array)[[2]] <- c("high", "low")
for (i in 1:length(pg.levels)){
  if (pg.data$Risk[1] == "high"){
    pg.add.array[i, "high", ] <- round(pg.data$Volume[1] * pg.data$Proportion * pg.levels[i], 0)
    pg.add.array[i, "low", ]  <- 0
  } else {
    pg.add.array[i, "high", ] <- 0
    pg.add.array[i, "low", ]  <- round(pg.data$Volume[1] * pg.data$Proportion * pg.levels[i], 0)
  }
}

## LOAD packages and functions
library(dplyr)
library(abind)
source("Profound-Function-PopInitialization.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-Function-Microsimulation.R")
source("Profound-DecisionTree.R")
source("Profound-DataInput.R")
source("Profound-Function-NxAvailAlgm.R")
source("Profound-CEA.R")
source("Profound-InputOutput-Setup.R")

## Initialize the study population - people who are at risk of opioid overdose
pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if(file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init_ppl  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
} else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init_ppl  <- pop.initiation(initials = initials, seed=seed)
  saveRDS(init_ppl, paste0("Inputs/InitialPopulation.rds"))
}

## Load calibrated parameters and seeds
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
sim.seed    <- sim.seed[1:100]     # ATTN: This is temporary, used to limit a smaller sample for analysis to save time

od.death.mx.last <- od.death.mx.totl <- matrix(0, nrow = length(sim.seed), ncol = 1+length(pg.levels))
scenario.name <- c("Status Quo", "100% increase", "500% increase", "1000% increase", "2000% increase", "5000% increase")
colnames(od.death.mx.last) <- colnames(od.death.mx.totl) <- scenario.name

for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  params <- sim.data.ls[[ss]] 
  params$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
  params$mor_Nx <- params$mor_bl * (1-0.9) # ATTN: This is temporary, used to manually force naloxone effect to be 90% (reducing probability of dying from an overdose)
  sim_sq          <- MicroSim(init_ppl, params = params, timesteps, v.state, d.c, PT.out = FALSE, Str = "SQ", seed = sim.seed[ss])        # run for status quo
  od.death.mx.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  od.death.mx.totl[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-59):timesteps, ])
  
  for (ll in 1:dim(pg.add.array)[1]){
    params$pg.add <- pg.add.array[ll, , ]
    sim_pg  <- MicroSim(init_ppl, params = params, timesteps, v.state, d.c, PT.out = FALSE, Str = "program", seed = sim.seed[ss]) # run for program scenario
    od.death.mx.last[ss, scenario.name[ll+1]] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    od.death.mx.totl[ss, scenario.name[ll+1]] <- sum(sim_pg$m.oddeath[(timesteps-59):timesteps, ])
  }
}

write.csv(od.death.mx.last, file = ("Outputs/Program/SA/SA_5Y_TotalODdeaths_last.csv"), row.names = F)
write.csv(od.death.mx.totl, file = ("Outputs/Program/SA/SA_5Y_TotalODdeaths_total.csv"), row.names = F)