

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

# LOAD packages and functions
library(openxlsx)
library(dplyr)
library(abind)
source("population.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_available.R")
source("cost_effectiveness.R")

# INPUT PARAMETERS
yr_start    <- 2016        # starting year of simulation 
yr_end     <- 2022        # end year of simulation (also the year for evaluation)
d.c         <- 0.03        # discounting of costs by 3%

source("io_setup.R")

## Initialize the study population - people who are at risk of opioid overdose
pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if(file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init_ppl  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
} else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init_ppl  <- pop.initiation(initials = initials, seed=seed)
  saveRDS(init_ppl, paste0("Inputs/InitialPopulation.rds"))
}

sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
sim.seed    <- sim.seed[1:50]

scenario.name <- c("Zero", "Status Quo", "Double", "Five times", "Ten times", "Saturation")
expand.level  <- c(0, 1, 2, 5, 10, 10000)
od.death.mx.last <- od.death.mx.wtns <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
colnames(od.death.mx.last) <- colnames(od.death.mx.wtns) <- scenario.name

for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  params.temp<- sim.data.ls[[ss]]
  params.temp$NxDataPharm$pe  <- 0
  params.temp$mor_Nx <- params.temp$mor_bl * (1-0.9)
  sim_sq          <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, Str = "SQ", seed = sim.seed[ss])        # run for status quo
  od.death.mx.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  od.death.mx.wtns[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-11):timesteps])
  
  exp.lv <- 0
  sim_pg  <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, Str = "expand", seed = sim.seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Zero"] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
  od.death.mx.wtns[ss, "Zero"] <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
  
  for (jj in 3:length(expand.level)){
    exp.lv  <- expand.level[jj]
    sim_pg  <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, Str = "expand", seed = sim.seed[ss]) # run for program scenario
    od.death.mx.last[ss, jj] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    od.death.mx.wtns[ss, jj] <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
  }
}

write.csv(od.death.mx.last, file = ("Ignore/SA/SA_Expand_3Y_TotalODdeaths.csv"), row.names = F)
write.csv(od.death.mx.wtns, file = ("Ignore/SA/SA_Expand_3Y_WitnessedDeaths.csv"), row.names = F)