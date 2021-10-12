#######################################################################################################
#####################            Sensitivity analysis - statewide program             #################
#######################################################################################################
# To assess the impact of statewide program (naloxone increase by the same multiplier in all regions)

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
yr.first     <- 2016        # starting year of simulation 
yr.last      <- 2022        # end year of simulation (also the year for evaluation)
d.c          <- 0.03        # discounting of costs by 3%
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

## LOAD packages and functions
library(openxlsx)
library(dplyr)
library(abind)
library(tictoc)
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
  init.pop  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
} else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init.pop  <- pop.initiation(initials = initials, seed=seed)
  saveRDS(init.pop, paste0("Inputs/InitialPopulation.rds"))
}

## Load calibrated parameters and seeds
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
sim.seed    <- sim.seed[1:50]

#define different statewide program expansion scenarios, including a 0 level and a saturation level
scenario.name <- c("Zero", "Status Quo", "10,000", "20,000", "50,000", "100,000", "150,000", "200,000", "250,000", "300,000", "350,000", "400,000", "450,000", "500,000", "600,000", "700,000", "800,000", "900,000", "1M")
sq.nlx     <- 8241
expand.nlx <- c(0, 8241, 10000, 20000, 50000, 100000, 150000, 200000, 250000, 300000, 350000, 400000, 450000, 500000, 600000, 700000, 800000, 900000, 1000000)
expand.level  <- expand.nlx / sq.nlx
od.death.mx.last <- od.death.mx.totl <- od.death.mx.wtns <- od.death.mx.wtns.tl <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
colnames(od.death.mx.last) <- colnames(od.death.mx.totl) <- colnames(od.death.mx.wtns) <- colnames(od.death.mx.wtns.tl) <- scenario.name

for (ss in 1:length(sim.seed)){
  # tic()
  print(paste0("Parameter set: ", ss))
  vparameters.temp<- sim.data.ls[[ss]]
  vparameters.temp$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
  vparameters.temp$mor_Nx <- vparameters.temp$mor_bl * (1-0.95) # ATTN: This is temporary, used to manually force naloxone effect to be 90% (reducing probability of dying from an overdose)
  sim_sq          <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "SQ", seed = sim.seed[ss])        # run for status quo
  od.death.mx.last[ss, "Status Quo"]    <- sum(sim_sq$m.oddeath[(n.t-11):n.t, ])
  od.death.mx.totl[ss, "Status Quo"]    <- sum(sim_sq$m.oddeath[(n.t-35):n.t, ])
  od.death.mx.wtns[ss, "Status Quo"]    <- sum(sim_sq$v.oddeath.w[(n.t-11):n.t])
  od.death.mx.wtns.tl[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(n.t-35):n.t])
  
  exp.lv <- 0
  sim_pg  <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "expand", seed = sim.seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Zero"]    <- sum(sim_pg$m.oddeath[(n.t-11):n.t, ])
  od.death.mx.totl[ss, "Zero"]    <- sum(sim_pg$m.oddeath[(n.t-35):n.t, ])
  od.death.mx.wtns[ss, "Zero"]    <- sum(sim_pg$v.oddeath.w[(n.t-11):n.t])
  od.death.mx.wtns.tl[ss, "Zero"] <- sum(sim_pg$v.oddeath.w[(n.t-35):n.t])
  
  for (jj in 3:length(expand.level)){
    exp.lv  <- expand.level[jj]
    sim_pg  <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "expand", seed = sim.seed[ss]) # run for program scenario
    od.death.mx.last[ss, jj]    <- sum(sim_pg$m.oddeath[(n.t-11):n.t, ])
    od.death.mx.totl[ss, jj]    <- sum(sim_pg$m.oddeath[(n.t-35):n.t, ])
    od.death.mx.wtns[ss, jj]    <- sum(sim_pg$v.oddeath.w[(n.t-11):n.t])
    od.death.mx.wtns.tl[ss, jj] <- sum(sim_pg$v.oddeath.w[(n.t-35):n.t])
  }
  # toc()
}

write.xlsx(od.death.mx.last, file = ("Ignore/Saturation/Saturation_3Y_TotalODdeaths_last.xlsx"), row.names = F)
write.xlsx(od.death.mx.totl, file = ("Ignore/Saturation/Saturation_3Y_TotalODdeaths_total.xlsx"), row.names = F)
write.xlsx(od.death.mx.wtns, file = ("Ignore/Saturation/Saturation_3Y_WitnessedDeaths_last.xlsx"), row.names = F)
write.xlsx(od.death.mx.wtns.tl, file = ("Ignore/Saturation/Saturation_3Y_WitnessedDeaths_total.xlsx"), row.names = F)