# TO_REVIEW does this need a separate file? Why not just use main with different inputs?

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

## Model setup parameters ##
seed <- 2021
sw.EMS.ODloc <- "overall" # Please choose from "overall" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "overall"

# Program data
library(openxlsx)
pg.data <- read.xlsx("Ignore/OEND_program.xlsx", sheet = "Project Weber")
pg.levels <- c(1, 5, 10, 20, 50)
pg.add.array <- array(0, dim = c(length(pg.levels), 2, dim(pg.data)[1]))
dimnames(pg.add.array)[[2]] <- c("high", "low")
for (i in 1:length(pg.levels)) {
  if (pg.data$Risk[1] == "high") {
    pg.add.array[i, "high", ] <- round(pg.data$Volume[1] * pg.data$Proportion * pg.levels[i], 0)
    pg.add.array[i, "low", ] <- 0
  } else {
    pg.add.array[i, "high", ] <- 0
    pg.add.array[i, "low", ] <- round(pg.data$Volume[1] * pg.data$Proportion * pg.levels[i], 0)
  }
}

# install.packages("rstudioapi")
# library(rstudioapi)
library(dplyr)
# library(tictoc)
library(abind)
source("population.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_available.R")
source("cost_effectiveness.R")

# INPUT PARAMETERS
yr_start <- 2016
yr_end <- 2024
pop.info <- c(
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
v.nlxused <- rep(0, times = timesteps) # count of naloxone kits used at each time step
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
pop.info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if (file.exists(paste0("Inputs/InitialPopulation.rds"))) {
  init_ppl <- readRDS(paste0("Inputs/InitialPopulation.rds"))
} else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))) {
  init_ppl <- initiate_ppl(initials = initials, seed = seed)
  saveRDS(init_ppl, paste0("Inputs/InitialPopulation.rds"))
}

simulation_data <- readRDS(file = paste0("calibration/CalibratedData.rds"))
simulation_seed <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
simulation_seed <- simulation_seed[1:100]

# sq.dh.mx  <- sq.nx.mx <- matrix(0, nrow = length(v.region), ncol = length(simulation_seed))
# pg.dh.ar  <- pg.nx.ar <- array(0, dim = c(dim(pg.add.array)[1], length(v.region), length(simulation_seed)))
# nlx.used.mx <- matrix(0, nrow = length(simulation_seed), ncol = 1+length(pg.levels))
od.death.mx.last <- od.death.mx.totl <- matrix(0, nrow = length(simulation_seed), ncol = 1 + length(pg.levels))
scenario.name <- c("Status Quo", "100% increase", "500% increase", "1000% increase", "2000% increase", "5000% increase")
# colnames(nlx.used.mx) <- scenario.name
colnames(od.death.mx.last) <- colnames(od.death.mx.totl) <- scenario.name

for (ss in 1:length(simulation_seed)) {
  print(paste0("Parameter set: ", ss))
  params.temp <- simulation_data[[ss]]
  params.temp$NxDataPharm$pe <- 0
  params.temp$mortality_nx <- params.temp$mor_bl * (1 - 0.9)
  sim_sq <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "SQ", seed = simulation_seed[ss]) # run for status quo
  # sq.dh.mx[ , ss] <- colSums(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  # sq.nx.mx[ , ss] <- colSums(sim_sq$n.nlx.OEND.str)
  # nlx.used.mx[ss, "Status Quo"] <- sum(sim_sq$v.nlxused[(timesteps-11):timesteps])
  od.death.mx.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps - 11):timesteps, ])
  od.death.mx.totl[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps - 59):timesteps, ])

  for (ll in 1:dim(pg.add.array)[1]) {
    params.temp$pg.add <- pg.add.array[ll, , ]
    sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "program", seed = simulation_seed[ss]) # run for program scenario
    # pg.dh.ar[ll, , ss] <- colSums(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    # pg.nx.ar[ll, , ss] <- colSums(sim_pg$n.nlx.OEND.str)
    # nlx.used.mx[ss, scenario.name[ll+1]] <- sum(sim_pg$v.nlxused[(timesteps-11):timesteps])
    od.death.mx.last[ss, scenario.name[ll + 1]] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])
    od.death.mx.totl[ss, scenario.name[ll + 1]] <- sum(sim_pg$m.oddeath[(timesteps - 59):timesteps, ])
  }
}

write.csv(od.death.mx.last, file = ("Ignore/SA/SA_5Y_TotalODdeaths_last.csv"), row.names = F)
write.csv(od.death.mx.totl, file = ("Ignore/SA/SA_5Y_TotalODdeaths_total.csv"), row.names = F)
