#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list=ls())
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

args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
# REVIEWED see analyze_calibration for decision on ov
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

# Program data
pg.data   <- read.xlsx("Ignore/OEND_program.xlsx", sheet = "Project Weber")
pg.levels <- c(1, 5, 10, 20, 50)
pg.add.array <- array(0, dim = c(length(pg.levels), 2, dim(pg.data)[1]))
# REVIEWED add.array = additional naloxone for high or low risk program
dimnames(pg.add.array)[[2]] <- c("high", "low")
for (i in 1:length(pg.levels)){
  # REVIEWED pg.levels = 1 means double... now is 2, 5, 10, 20, 50; maybe use risk_level instead of risk (risk sounds like a prob)
  if (pg.data$Risk[1] == "high"){
    pg.add.array[i, "high", ] <- round(pg.data$Volume[1] * pg.data$Proportion * pg.levels[i], 0)
    pg.add.array[i, "low", ]  <- 0
  } else {
    pg.add.array[i, "high", ] <- 0
    pg.add.array[i, "low", ]  <- round(pg.data$Volume[1] * pg.data$Proportion * pg.levels[i], 0)
  }
}


# INPUT PARAMETERS
yr_start    <- 2016
yr_end     <- 2020
d.c         <- 0.03   

source("io_setup.R")

## Initialize the study population - people who are at risk of opioid overdose
ppl_info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if(file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init_ppl  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
} else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))){
  init_ppl  <- ppl.initiation(initials = initials, seed=seed)
  saveRDS(init_ppl, paste0("Inputs/InitialPopulation.rds"))
}

# REVIEWED ls = list
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
sim.seed    <- sim.seed[1:100]

# REVIEWED sq = status quo; pg = program
sq.dh.mx  <- sq.nx.mx <- matrix(0, nrow = length(v.region), ncol = length(sim.seed))
pg.dh.ar  <- pg.nx.ar <- array(0, dim = c(dim(pg.add.array)[1], length(v.region), length(sim.seed)))
nlx.used.mx <- matrix(0, nrow = length(sim.seed), ncol = 1+length(pg.levels))
od.death.mx <- matrix(0, nrow = length(sim.seed), ncol = 1+length(pg.levels))
scenario.name <- c("Status Quo", "100% increase", "500% increase", "1000% increase", "2000% increase", "5000% increase")
colnames(nlx.used.mx) <- scenario.name
colnames(od.death.mx) <- scenario.name

for (ss in 1:length(sim.seed)){
  print(paste0("Parameter set: ", ss))
  params.temp<- sim.data.ls[[ss]]
  params.temp$NxDataPharm$pe  <- 0
  params.temp$mor_Nx <- params.temp$mor_bl * (1-0.9)
  sim_sq          <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "SQ", seed = sim.seed[ss])        # run for status quo
  sq.dh.mx[ , ss] <- colSums(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  sq.nx.mx[ , ss] <- colSums(sim_sq$n.nlx.OEND.str)
  nlx.used.mx[ss, "Status Quo"] <- sum(sim_sq$v.nlxused[(timesteps-11):timesteps])
  od.death.mx[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  
  for (ll in 1:dim(pg.add.array)[1]){
    params.temp$pg.add <- pg.add.array[ll, , ]
    sim_pg  <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "program", seed = sim.seed[ss]) # run for program scenario
    pg.dh.ar[ll, , ss] <- colSums(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    pg.nx.ar[ll, , ss] <- colSums(sim_pg$n.nlx.OEND.str)
    nlx.used.mx[ss, scenario.name[ll+1]] <- sum(sim_pg$v.nlxused[(timesteps-11):timesteps])
    od.death.mx[ss, scenario.name[ll+1]] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
  }
}

pop.region  <- colSums(Demographic[ , -c(1:3)])

preliminary.NoDeaths <- data.frame(matrix(nrow = num_regions * (1+dim(pg.add.array)[1]), ncol = 5))
x <- c("location", "scenario", "mean", "upper", "lower")
colnames(preliminary.NoDeaths) <- x
preliminary.NoDeaths$location <- rep(v.region, 1+dim(pg.add.array)[1])

preliminary.NoDeaths$scenario <- rep(scenario.name, each  = length(v.region))

preliminary.RateNlx <- preliminary.NoNlx <- preliminary.RateDeaths <- preliminary.NoDeaths
#Number of deaths
preliminary.NoDeaths$mean[preliminary.NoDeaths$scenario == "Status Quo"]  <- apply(sq.dh.mx, 1, mean)
preliminary.NoDeaths$upper[preliminary.NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, quantile, probs = 0.975)
preliminary.NoDeaths$lower[preliminary.NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, quantile, probs = 0.025)

for (sc in 2:length(scenario.name)){
  preliminary.NoDeaths$mean[preliminary.NoDeaths$scenario == scenario.name[sc]]  <- apply(pg.dh.ar[sc-1,,], 1, mean)
  preliminary.NoDeaths$upper[preliminary.NoDeaths$scenario == scenario.name[sc]] <- apply(pg.dh.ar[sc-1,,], 1, quantile, probs = 0.975)
  preliminary.NoDeaths$lower[preliminary.NoDeaths$scenario == scenario.name[sc]] <- apply(pg.dh.ar[sc-1,,], 1, quantile, probs = 0.025)
}

#Rate of deaths
preliminary.RateDeaths[ , c("mean", "upper", "lower")] <- preliminary.NoDeaths[ , c("mean", "upper", "lower")] / pop.region * 100000

#Number of Naloxone kits
preliminary.NoNlx$mean[preliminary.NoNlx$scenario == "Status Quo"]  <- apply(sq.nx.mx, 1, mean)
preliminary.NoNlx$upper[preliminary.NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, quantile, probs = 0.975)
preliminary.NoNlx$lower[preliminary.NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, quantile, probs = 0.025)

for (sc in 2:length(scenario.name)){
  preliminary.NoNlx$mean[preliminary.NoNlx$scenario == scenario.name[sc]]  <- apply(pg.nx.ar[sc-1,,], 1, mean)
  preliminary.NoNlx$upper[preliminary.NoNlx$scenario == scenario.name[sc]] <- apply(pg.nx.ar[sc-1,,], 1, quantile, probs = 0.975)
  preliminary.NoNlx$lower[preliminary.NoNlx$scenario == scenario.name[sc]] <- apply(pg.nx.ar[sc-1,,], 1, quantile, probs = 0.025)
}

#Rate of Naloxone kits
preliminary.RateNlx[ , c("mean", "upper", "lower")] <- preliminary.NoNlx[ , c("mean", "upper", "lower")] / pop.region * 100000

write.csv(preliminary.NoDeaths, file = ("Ignore/preliminary.Number.Deaths.csv"), row.names = F)
write.csv(preliminary.RateDeaths, file = ("Ignore/preliminary.Rate.Deaths.csv"), row.names = F)
write.csv(preliminary.NoNlx, file = ("Ignore/preliminary.Number.Naloxone.csv"), row.names = F)
write.csv(preliminary.RateNlx, file = ("Ignore/preliminary.Rate.Naloxone.csv"), row.names = F)
write.csv(nlx.used.mx, file = ("Ignore/preliminary.NaloxoneUsed.csv"), row.names = F)
write.csv(od.death.mx, file = ("Ignore/preliminary.TotalODdeaths.csv"), row.names = F)
