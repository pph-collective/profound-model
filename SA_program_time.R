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

# Program data
library(openxlsx)
pg.data   <- read.xlsx("Ignore/OEND_program.xlsx", sheet = "Project Weber") #load data for the program, including high/low-risk and last year's numbers
pg.levels <- c(2, 5, 10, 20, 50)   #define program expansion levels
scenario.name <- c("Status Quo", "Double", "5 times", "10 times", "20 times", "50 times")
pg.add.array <- array(0, dim = c(length(pg.levels), 2, dim(pg.data)[1]))     # array for the added naloxone kits (only acount for the additional kits) due to program expansion, stratified by high/low-risk program and regions
dimnames(pg.add.array)[[2]] <- c("high", "low")
for (i in 1:length(pg.levels)){
  if (pg.data$Risk[1] == "high"){
    pg.add.array[i, "high", ] <- round(pg.data$Volume[1] * pg.data$Proportion * (pg.levels[i]-1), 0)
    pg.add.array[i, "low", ]  <- 0
  } else {
    pg.add.array[i, "high", ] <- 0
    pg.add.array[i, "low", ]  <- round(pg.data$Volume[1] * pg.data$Proportion * (pg.levels[i]-1), 0)
  }
}

## LOAD packages and functions
library(openxlsx)
library(dplyr)
library(abind)
source("Profound-Function-PopInitialization.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-Function-Microsimulation(Trial).R")
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

m.oddeath.preb <- m.oddeath.il.lr <- m.oddeath.il.hr <- m.nlx.mn <- rep(0, times = n.t)
od.death.mx.first <- od.death.mx.total <- od.death.mx.last <- matrix(0, nrow = length(sim.seed), ncol = 1+length(pg.levels))
colnames(od.death.mx.first) <- colnames(od.death.mx.total) <- colnames(od.death.mx.last) <- scenario.name
v.nlx.mn <- matrix(0, nrow= n.t, ncol = length(scenario.name))
colnames(v.nlx.mn) <- scenario.name
for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  vparameters.temp<- sim.data.ls[[ss]]
  vparameters.temp$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
  vparameters.temp$mor_Nx <- vparameters.temp$mor_bl * (1-0.9) # ATTN: This is temporary, used to manually force naloxone effect to be 90% (reducing probability of dying from an overdose)
  sim_sq          <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "SQ", seed = sim.seed[ss])        # run for status quo
  v.nlx.mn[ , "Status Quo"] <- sim_sq$m.nlx.mn
  od.death.mx.first[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[49:60, ])
  od.death.mx.last[ss, "Status Quo"]  <- sum(sim_sq$m.oddeath[(n.t-11):n.t, ])
  od.death.mx.total[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[49:n.t, ])
  
  for (ll in 1:dim(pg.add.array)[1]){
    vparameters.temp$pg.add <- pg.add.array[ll, , ]
    sim_pg  <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "program", seed = sim.seed[ss]) # run for program scenario
    v.nlx.mn[ , ll+1] <- sim_pg$m.nlx.mn
    od.death.mx.first[ss, ll+1] <- sum(sim_pg$m.oddeath[49:60, ])
    od.death.mx.last[ss, ll+1]  <- sum(sim_pg$m.oddeath[(n.t-11):n.t, ])
    od.death.mx.total[ss, ll+1] <- sum(sim_pg$m.oddeath[49:n.t, ])
  }
}
write.xlsx(od.death.mx.first, file="Ignore/SA/time/Program3Yfirst.xlsx",
           col.names = T, row.names = F)
write.xlsx(od.death.mx.total, file="Ignore/SA/time/Program3Ytotal.xlsx",
           col.names = T, row.names = F)
write.xlsx(od.death.mx.last, file="Ignore/SA/time/Program3Ylast.xlsx",
           col.names = T, row.names = F)
write.xlsx(v.nlx.mn, file="Ignore/SA/time/NaloxoneMonthly.xlsx",
           col.names = T, row.names = F)
