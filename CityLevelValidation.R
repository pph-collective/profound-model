#################################################################################################
###################                Model validation - city-level         ########################
#################################################################################################

###########################################################################################################
####    Additional validation using calibrated parameters against targets at jurisdiction level        ####
####    1 target  -  annual # overdose deaths in each jurisdiction:       od.deathYR (YR for year)     ####
####    Need to rerun model over selected sets of calibrated parameters                                ####
####    Rely on many inputs and functions loaded in Calibration_ResultAnalysis.R                       ####
###########################################################################################################

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
yr.first     <- 2016   #first year of evaluation
yr.last      <- 2020   #last year of evaluation
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

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
  init.pop  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
} else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init.pop  <- pop.initiation(initials = initials, seed=seed)
  saveRDS(init.pop, paste0("Inputs/InitialPopulation.rds"))
}

## Load calibrated parameters and seeds
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
sim.seed    <- sim.seed[1:50]

# Initiate result matrices
ODdeaths16 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
ODdeaths17 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
ODdeaths18 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
ODdeaths19 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
row.names(ODdeaths16) <- row.names(ODdeaths17) <- row.names(ODdeaths18) <- row.names(ODdeaths19) <- sim.data.ls[[1]]$v.rgn

for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  vparameters.temp<- sim.data.ls[[ss]]
  sim_sq    <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "SQ", seed = sim.seed[ss])        # run for status quo
  ODdeaths16[ , ss] <- colSums(sim_sq$m.oddeath[1:12, ])
  ODdeaths17[ , ss] <- colSums(sim_sq$m.oddeath[13:24, ])
  ODdeaths18[ , ss] <- colSums(sim_sq$m.oddeath[25:36, ])
  ODdeaths19[ , ss] <- colSums(sim_sq$m.oddeath[37:48, ])
}

write.xlsx(ODdeaths16, file="Outputs/CityLevelValidatin.xlsx", sheetName = "2016",
           col.names = F, row.names = T)
write.xlsx(ODdeaths17, file="Outputs/CityLevelValidatin.xlsx", sheetName = "2017", append=TRUE,
           col.names = F, row.names = T)
write.xlsx(ODdeaths18, file="Outputs/CityLevelValidatin.xlsx", sheetName = "2018", append=TRUE,
           col.names = F, row.names = T)
write.xlsx(ODdeaths19, file="Outputs/CityLevelValidatin.xlsx", sheetName = "2019", append=TRUE,
           col.names = F, row.names = T)