###########################################################################################################
#####################            Program evaluation - individual program             ######################
###########################################################################################################
# To evaluate the impact of specific OEND program at a few hypothetical expansion levels

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
yr.first     <- 2016   #first year of evaluation
yr.last      <- 2020   #last year of evaluation
d.c          <- 0.03   #discounting rate (for costs)
seed         <- 2021   #initial seed for stochastic sampling
sw.EMS.ODloc <- "ov"   #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

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
sim.seed    <- sim.seed[1:100]    # ATTN: This is temporary, used to limit a smaller sample for analysis to save time

sq.dh.mx  <- sq.nx.mx <- matrix(0, nrow = length(v.rgn), ncol = length(sim.seed))                  #status quo, od death matrix, rows for regions and columns for simulations
pg.dh.ar  <- pg.nx.ar <- array(0, dim = c(dim(pg.add.array)[1], length(v.rgn), length(sim.seed)))  #program scenario, od death array, dimensions are program level, region, simulations
nlx.used.mx <- matrix(0, nrow = length(sim.seed), ncol = 1+length(pg.levels))                      #matrix for annual # naloxone kits distributed in the last evaluation year, rows for simulations and columns for scenarios (including status quo)            
od.death.mx <- matrix(0, nrow = length(sim.seed), ncol = 1+length(pg.levels))                      #matrix for state-level od deaths, rows for simulations and columns for scenarios (including status quo)            
colnames(nlx.used.mx) <- scenario.name
colnames(od.death.mx) <- scenario.name

for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  vparameters.temp<- sim.data.ls[[ss]]
  vparameters.temp$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
  vparameters.temp$mor_Nx <- vparameters.temp$mor_bl * (1-0.9) # ATTN: This is temporary, used to manually force naloxone effect to be 90% (reducing probability of dying from an overdose)
  sim_sq          <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "SQ", seed = sim.seed[ss])        # run for status quo
  sq.dh.mx[ , ss] <- colSums(sim_sq$m.oddeath[(n.t-11):n.t, ])
  sq.nx.mx[ , ss] <- colSums(sim_sq$n.nlx.OEND.str)
  nlx.used.mx[ss, "Status Quo"] <- sum(sim_sq$v.nlxused[(n.t-11):n.t])
  od.death.mx[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(n.t-11):n.t, ])
  
  for (ll in 1:dim(pg.add.array)[1]){
    vparameters.temp$pg.add <- pg.add.array[ll, , ]
    sim_pg  <- MicroSim(init.pop, vparameters = vparameters.temp, n.t, v.state, d.c, PT.out = FALSE, Str = "program", seed = sim.seed[ss]) # run for program scenario
    pg.dh.ar[ll, , ss] <- colSums(sim_pg$m.oddeath[(n.t-11):n.t, ])
    pg.nx.ar[ll, , ss] <- colSums(sim_pg$n.nlx.OEND.str)
    nlx.used.mx[ss, scenario.name[ll+1]] <- sum(sim_pg$v.nlxused[(n.t-11):n.t])
    od.death.mx[ss, scenario.name[ll+1]] <- sum(sim_pg$m.oddeath[(n.t-11):n.t, ])
  }
}

pop.rgn  <- colSums(Demographic[ , -c(1:3)])   #population in each region (for calculation of rate)

NoDeaths <- data.frame(matrix(nrow = n.rgn * (1+dim(pg.add.array)[1]), ncol = 5))
x <- c("location", "scenario", "mean", "upper", "lower")
colnames(NoDeaths) <- x
NoDeaths$location <- rep(v.rgn, 1+dim(pg.add.array)[1])

NoDeaths$scenario <- rep(scenario.name, each  = length(v.rgn))

RateNlx <- NoNlx <- RateDeaths <- NoDeaths
#Number of deaths
NoDeaths$mean[NoDeaths$scenario == "Status Quo"]  <- apply(sq.dh.mx, 1, mean)
NoDeaths$upper[NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, quantile, probs = 0.975)
NoDeaths$lower[NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, quantile, probs = 0.025)

for (sc in 2:length(scenario.name)){
  NoDeaths$mean[NoDeaths$scenario == scenario.name[sc]]  <- apply(pg.dh.ar[sc-1,,], 1, mean)
  NoDeaths$upper[NoDeaths$scenario == scenario.name[sc]] <- apply(pg.dh.ar[sc-1,,], 1, quantile, probs = 0.975)
  NoDeaths$lower[NoDeaths$scenario == scenario.name[sc]] <- apply(pg.dh.ar[sc-1,,], 1, quantile, probs = 0.025)
}

#Rate of deaths
RateDeaths[ , c("mean", "upper", "lower")] <- NoDeaths[ , c("mean", "upper", "lower")] / pop.rgn * 100000

#Number of Naloxone kits
NoNlx$mean[NoNlx$scenario == "Status Quo"]  <- apply(sq.nx.mx, 1, mean)
NoNlx$upper[NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, quantile, probs = 0.975)
NoNlx$lower[NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, quantile, probs = 0.025)

for (sc in 2:length(scenario.name)){
  NoNlx$mean[NoNlx$scenario == scenario.name[sc]]  <- apply(pg.nx.ar[sc-1,,], 1, mean)
  NoNlx$upper[NoNlx$scenario == scenario.name[sc]] <- apply(pg.nx.ar[sc-1,,], 1, quantile, probs = 0.975)
  NoNlx$lower[NoNlx$scenario == scenario.name[sc]] <- apply(pg.nx.ar[sc-1,,], 1, quantile, probs = 0.025)
}

#Rate of Naloxone kits
RateNlx[ , c("mean", "upper", "lower")] <- NoNlx[ , c("mean", "upper", "lower")] / pop.rgn * 100000

write.csv(NoDeaths, file = ("Outputs/Program/Number.Deaths.csv"), row.names = F)
write.csv(RateDeaths, file = ("Outputs/Program//Rate.Deaths.csv"), row.names = F)
write.csv(NoNlx, file = ("Outputs/Program//Number.Naloxone.csv"), row.names = F)
write.csv(RateNlx, file = ("Outputs/Program//Rate.Naloxone.csv"), row.names = F)
write.csv(nlx.used.mx, file = ("Outputs/Program//NaloxoneUsed.csv"), row.names = F)
write.csv(od.death.mx, file = ("Outputs/Program//TotalODdeaths.csv"), row.names = F)