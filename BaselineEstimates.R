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
yr_start     <- 2016        # starting year of simulation 
yr_end       <- 2020        # end year of simulation (also the year for evaluation)
discount.rate<- 0.03        # discounting of costs by 3%
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

## LOAD packages and functions
library(xlsx)
library(dplyr)
library(abind)
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")
source("io_setup.R")

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
sim.seed    <- sim.seed[1:50]

m.oddeath.preb <- m.oddeath.il.lr <- m.oddeath.il.hr <- m.nlx.mn <- rep(0, times = timesteps)
od.death.mx    <- matrix(0, nrow= length(sim.seed), ncol = 4)
colnames(od.death.mx) <- c("preb","il.lr", "il.hr", "stim")
pop.mx    <- matrix(0, nrow= length(sim.seed), ncol = 5)
colnames(pop.mx) <- c("preb","il.lr", "il.hr", "stim", "inact")
v.nlx.mn <- matrix(0, nrow= length(sim.seed), ncol = timesteps)
for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  params.temp<- sim.data.ls[[ss]]
  params.temp$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
  # vparameters.temp$mor_Nx <- vparameters.temp$mor_bl * (1-0.9) # ATTN: This is temporary, used to manually force naloxone effect to be 90% (reducing probability of dying from an overdose)
  sim_sq          <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, discount.rate, PT.out = TRUE, strategy = "SQ", seed = sim.seed[ss])        # run for status quo
  od.death.mx[ss, "preb"]  <- sum(sim_sq$m.oddeath.preb[37:48])
  od.death.mx[ss, "il.lr"] <- sum(sim_sq$m.oddeath.il.lr[37:48])
  od.death.mx[ss, "il.hr"] <- sum(sim_sq$m.oddeath.il.hr[37:48])
  od.death.mx[ss, "stim"]  <- sum(sim_sq$m.oddeath.st[37:48])
  
  pop.mx[ss, "preb"] <- nrow(subset(sim_sq$pop.trace[[42]], curr.state == "relap" & OU.state == "preb")) + nrow(subset(sim_sq$pop.trace[[42]], curr.state == "preb"))
  pop.mx[ss, "il.lr"] <- nrow(subset(sim_sq$pop.trace[[42]], curr.state == "relap" & OU.state == "il.lr")) + nrow(subset(sim_sq$pop.trace[[42]], curr.state == "il.lr"))
  pop.mx[ss, "il.hr"] <- nrow(subset(sim_sq$pop.trace[[42]], curr.state == "relap" & OU.state == "il.hr")) + nrow(subset(sim_sq$pop.trace[[42]], curr.state == "il.hr"))
  pop.mx[ss, "stim"] <- nrow(subset(sim_sq$pop.trace[[42]], curr.state == "NODU"))
  pop.mx[ss, "inact"] <- nrow(subset(sim_sq$pop.trace[[42]], curr.state == "inact"))
  v.nlx.mn[ss, ] <- sim_sq$m.nlx.mn
}
write.xlsx(od.death.mx, file="Ignore/BaselineEstimates.xlsx",
           col.names = T, row.names = F)
write.xlsx(pop.mx, file="Ignore/BaselinePopEstimates.xlsx",
           col.names = T, row.names = F)
write.xlsx(v.nlx.mn, file="Ignore/NaloxoneMonthly.xlsx",
           col.names = F, row.names = F)
