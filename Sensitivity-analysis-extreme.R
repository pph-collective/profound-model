#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list = ls())
args <- commandArgs(trailingOnly = TRUE)

## Model setup parameters ##
seed <- 2021
sw.EMS.ODloc <- "overall" # Please choose from "overall" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "overall"

# # Program data
library(openxlsx)

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
yr_end <- 2022
pop.info <- c(
  "sex", "race", "age", "residence", "curr.state",
  "OU.state", "init.age", "init.state", "ever.od", "fx"
) # information for each model individual
agent_states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead") # vector for state names
v.oustate <- c("preb", "il.lr", "il.hr") # vector for active opioid use state names
num_states <- length(agent_states) # number of states
num_years <- yr_end - yr.first + 1
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
simulation_seed <- simulation_seed[1:10]

# sq.dh.mx  <- sq.nx.mx <- matrix(0, nrow = length(v.region), ncol = length(simulation_seed))
# pg.dh.ar  <- pg.nx.ar <- array(0, dim = c(dim(pg.add.array)[1], length(v.region), length(simulation_seed)))
# nlx.used.mx <- matrix(0, nrow = length(simulation_seed), ncol = 1+length(pg.levels))
scenario.name <- c("Zero", "Status Quo", "Double", "Five times", "Ten times", "Saturation")
od.death.mx.last <- od.death.mx.totl <- matrix(0, nrow = length(simulation_seed), ncol = length(scenario.name))

colnames(od.death.mx.last) <- colnames(od.death.mx.totl) <- scenario.name

for (ss in 1:length(simulation_seed)) {
  print(paste0("Parameter set: ", ss))
  params.temp <- simulation_data[[ss]]
  params.temp$NxDataPharm$pe <- 0
  params.temp$mortality_nx <- params.temp$mor_bl * (1 - 0.9)
  sim_sq <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "SQ", seed = simulation_seed[ss]) # run for status quo

  od.death.mx.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps - 11):timesteps, ])

  exp.lv <- 0
  sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "expand", seed = simulation_seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Zero"] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])

  exp.lv <- 2
  sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, Str = "expand", seed = simulation_seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Double"] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])

  exp.lv <- 5
  sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, Str = "expand", seed = simulation_seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Five times"] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])

  exp.lv <- 10
  sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, Str = "expand", seed = simulation_seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Ten times"] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])

  exp.lv <- 10000
  sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "expand", seed = simulation_seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Saturation"] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])
}

write.csv(od.death.mx.last, file = ("Ignore/SA/SA_Expand_3Y_TotalODdeaths_last.csv"), row.names = F)
