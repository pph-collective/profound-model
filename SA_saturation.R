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
yr_start <- 2016        # starting year of simulation 
yr_end <- 2022        # end year of simulation (also the year for evaluation)
discount.rate <- 0.03        # discounting of costs by 3%
# seed <- 2021

## LOAD packages and functions
library("argparser")
library(openxlsx)
library(dplyr)
library(abind)
library(tictoc)
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")

args <- arg_parser("arguments")
args <- add_argument(args, "--seed", help = "seed for random numbers", default = 2021)
args <- add_argument(args, "--regional", help = "flag to run regional model", flag = TRUE)
args <- add_argument(args, "--outfile", help = "file to store outputs", default = "OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--ppl", help = "file with initial ppl info", default = "Inputs/init_pop.rds")
argv <- parse_args(args)
seed <- as.integer(argv$seed)
init_ppl.file <- argv$ppl
source("io_setup.R")

# ## Initialize the study population - people who are at risk of opioid overdose
# pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
# if(file.exists(paste0("Inputs/InitialPopulation.rds"))){
#   init.pop  <- readRDS(paste0("Inputs/InitialPopulation.rds"))
# } else if (!file.exists(paste0("Inputs/InitialPopulation.rds"))){
#   init.pop  <- pop.initiation(initials = initials, seed=seed)
#   saveRDS(init.pop, paste0("Inputs/InitialPopulation.rds"))
# }

## Load calibrated parameters and seeds
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
sim.seed    <- sim.seed[1:100]

#define different statewide program expansion scenarios, including a 0 level and a saturation level
# scenario.name <- c("Zero", "Status Quo", "10,000", "20,000", "30,000", "40,000", "50,000", "60,000", "70,000", "80,000", "90,000", "100,000")
sq.nlx     <- 8241
# expand.nlx <- c(0, 8241, 10000, 20000, 30000, 40000, 50000, 60000, 70000, 80000, 90000, 100000)
scenario.name <- c("Zero", "Status Quo", "200,000", "400,000","600,000", "800,000", "1,000,000")
expand.nlx <- c(0, 8241, 200000, 400000, 600000, 800000, 1000000)
expand.level  <- expand.nlx / sq.nlx
od.death.mx.last <- od.death.mx.totl <- od.death.mx.wtns <- od.death.mx.wtns.tl <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
colnames(od.death.mx.last) <- colnames(od.death.mx.totl) <- colnames(od.death.mx.wtns) <- colnames(od.death.mx.wtns.tl) <- scenario.name

for (ss in 1:length(sim.seed)){
  # tic()
  print(paste0("Parameter set: ", ss))
  vparameters.temp<- sim.data.ls[[ss]]
  # vparameters.temp$NxDataPharm$pe  <- 0 # ATTN: This is temporary, used to manually remove pharamacy naloxone
  # vparameters.temp$mor_Nx <- vparameters.temp$mor_bl * (1-0.95) # ATTN: This is temporary, used to manually force naloxone effect to be 90% (reducing probability of dying from an overdose)
  sim_sq <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = sim.seed[ss])        # run for status quo
  od.death.mx.last[ss, "Status Quo"]    <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  od.death.mx.totl[ss, "Status Quo"]    <- sum(sim_sq$m.oddeath[(timesteps-35):timesteps, ])
  od.death.mx.wtns[ss, "Status Quo"]    <- sum(sim_sq$v.oddeath.w[(timesteps-11):timesteps])
  od.death.mx.wtns.tl[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-35):timesteps])
  
  exp.lv <- 0
  sim_pg <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "expand", seed = sim.seed[ss]) # run for program scenario
  od.death.mx.last[ss, "Zero"]    <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
  od.death.mx.totl[ss, "Zero"]    <- sum(sim_pg$m.oddeath[(timesteps-35):timesteps, ])
  od.death.mx.wtns[ss, "Zero"]    <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
  od.death.mx.wtns.tl[ss, "Zero"] <- sum(sim_pg$v.oddeath.w[(timesteps-35):timesteps])
  
  for (jj in 3:length(expand.level)){
    exp.lv  <- expand.level[jj]
    sim_pg <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "expand", seed = sim.seed[ss]) # run for program scenario
    od.death.mx.last[ss, jj]    <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    od.death.mx.totl[ss, jj]    <- sum(sim_pg$m.oddeath[(timesteps-35):timesteps, ])
    od.death.mx.wtns[ss, jj]    <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
    od.death.mx.wtns.tl[ss, jj] <- sum(sim_pg$v.oddeath.w[(timesteps-35):timesteps])
  }
  # toc()
}

write.xlsx(od.death.mx.last, file = ("Ignore/Saturation/Saturation_3Y_TotalODdeaths_last2.xlsx"), row.names = F)
write.xlsx(od.death.mx.totl, file = ("Ignore/Saturation/Saturation_3Y_TotalODdeaths_total2.xlsx"), row.names = F)
write.xlsx(od.death.mx.wtns, file = ("Ignore/Saturation/Saturation_3Y_WitnessedDeaths_last2.xlsx"), row.names = F)
write.xlsx(od.death.mx.wtns.tl, file = ("Ignore/Saturation/Saturation_3Y_WitnessedDeaths_total2.xlsx"), row.names = F)