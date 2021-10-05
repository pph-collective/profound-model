#######################################################################################################
#####################            Sensitivity analysis - statewide program             #################
#######################################################################################################
# To assess the impact of statewide program (naloxone increase by the same multiplier in all regions)

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list = ls())

## Model setup parameters ##
seed <- 2021
sw.EMS.ODloc <- "overall" # Please choose from "overall" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "overall"

## LOAD packages and functions
source("main.R")
library(commandArgs)

# # INPUT PARAMETERS
# args <- commandArgs(trailingOnly = TRUE)
# year_start <- 2016 # starting year of simulation
# yr_end <- 2022 # end year of simulation (also the year for evaluation)
# d.c <- 0.03 # discounting of costs by 3%

# source("io_setup.R")

## Initialize the study population - people who are at risk of opioid overdose
pop.info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever_od", "fx")
if (file.exists(paste0("Inputs/init_pop.rds"))) {
  init_ppl <- readRDS(paste0("Inputs/init_pop.rds"))
} else if (!file.exists(paste0("Inputs/init_pop.rds"))) {
  init_ppl <- pop.initiation(initials = initials, seed = seed)
  saveRDS(init_ppl, paste0("Inputs/init_pop.rds"))
}

simulation_data <- readRDS(file = paste0("calibration/CalibratedData.rds"))
simulation_seed <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
simulation_seed <- simulation_seed[1:50]

# define different statewide program expansion scenarios, including a 0 level and a saturation level
scenario.name <- c("Zero", "Status Quo", "Double", "Five times", "Ten times", "Saturation")
expand.level <- c(0, 1, 2, 5, 10, 10000)
od.death.mx.last <- od.death.mx.wtns <- matrix(0, nrow = length(simulation_seed), ncol = length(scenario.name))
colnames(od.death.mx.last) <- colnames(od.death.mx.wtns) <- scenario.name

for (ss in 1:length(simulation_seed)) {
  # is this all in main??
  print(paste0("Parameter set: ", ss))
  params.temp <- simulation_data[[ss]]
  # ATTN: These two lines are temporary, to remove pharm nx and make nx effect 90%
  params.temp$nlx_data_pharm$pe <- 0
  params.temp$mortality_nx <- params.temp$mor_bl * (1 - 0.9)
  sim_sq <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "SQ", seed = simulation_seed[ss]) # run for status quo
  od.death.mx.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps - 11):timesteps, ])
  od.death.mx.wtns[ss, "Status Quo"] <- sum(sim_sq$death_wtns[(timesteps - 11):timesteps])

  exp.lv <- 0
  sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "expand", seed = simulation_seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Zero"] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])
  od.death.mx.wtns[ss, "Zero"] <- sum(sim_pg$death_wtns[(timesteps - 11):timesteps])

  for (jj in 3:length(expand.level)) {
    exp.lv <- expand.level[jj]
    sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "expand", seed = simulation_seed[ss]) # run for program scenario
    od.death.mx.last[ss, jj] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])
    od.death.mx.wtns[ss, jj] <- sum(sim_pg$death_wtns[(timesteps - 11):timesteps])
  }
}

write.csv(od.death.mx.last, file = ("Ignore/SA/SA_Expand_3Y_TotalODdeaths.csv"), row.names = F)
write.csv(od.death.mx.wtns, file = ("Ignore/SA/SA_Expand_3Y_WitnessedDeaths.csv"), row.names = F)
