###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2020 #########################
###############################################################################################
# Main module for the microsimulation of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: May 06, 2020
# Last update: April 08, 2021
#
###############################################################################################
#######################           Microsimulation           ###################################
###############################################################################################

###############################################################################################
####    Individual-based microsimulation                                                   ####
####    6 health states: prescribed, illicit (L/H), inactive, non-opioid, relapsed, death  ####
####    1 health event:  Overdose                                                          ####
####    Attributes:      state, age, sex, fentanyl, overdosed, pre.state,                  ####
####    Built to inform Naloxone distribution strategies to prevent overdsoe death         ####
###############################################################################################

#############################################################################
# 1. SET directpry and workspace
#############################################################################
rm(list=ls())
library("argparser")
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(tictoc)
source("population.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_available.R")
source("cost_effectiveness.R")

args = arg_parser("arguments")
args <- add_argument(args, "--seed", help="seed for random number", default=2021)
args <- add_argument(args, "--regional", help="flag to run regional model", flag=TRUE)
args <- add_argument(args, "--outfile", help="file to store outputs", default="OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--ppl", help="file with initial ppl info", default="Inputs/InitialPopulation.rds")
argv <- parse_args(args)
seed <- as.integer(argv$seed)

init.pop.file <- argv$ppl
## Model setup parameters ##
sw.EMS.ODloc <- "ov"
out.file <- argv$outfile
if (isTRUE(argv$regional)){
  sw.EMS.ODloc <- "sp"
}



# INPUT PARAMETERS
yr_start    <- 2016
yr_end     <- 2020
ppl_info    <- c("sex", "race", "age", "residence", "curr.state",
                 "OU.state", "init.age", "init.state", "ever.od", "fx")            # information for each model individual
agent_states     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
num_states     <- length(agent_states)                                                     # number of states
num_years        <- yr_end-yr_start+1
timesteps         <- 12 * num_years                                                           # number of time cycles (in month)
n.region       <- length(v.region)                                                       # number of regions

# OUTPUT matrices and vectors
# REVIEWED now in separate script
v.od        <- rep(0, times = timesteps)                                                 # count of overdose events at each time step
v.oddeath   <- rep(0, times = timesteps)                                                 # count of overdose deaths at each time step
m.oddeath   <- matrix(0, nrow = timesteps, ncol = n.region)
colnames(m.oddeath) <- v.region
v.odpriv    <- rep(0, times = timesteps)                                                 # count of overdose events occurred at private setting at each time step
v.odpubl    <- rep(0, times = timesteps)                                                 # count of overdose events occurred at public setting at each time step
v.deathpriv <- rep(0, times = timesteps)                                                 # count of overdose deaths occurred at private setting at each time step
v.deathpubl <- rep(0, times = timesteps)                                                 # count of overdose deaths occurred at public setting at each time step
v.nlxused   <- rep(0, times = timesteps)                                                 # count of naloxone kits used at each time step
v.str       <- c("SQ", "expand")                                                   # store the strategy names
d.c         <- 0.03                                                                # discounting of costs by 3%
cost.item   <- c("TotalCost", "NxCost")
cost.matrix <- matrix(0, nrow=timesteps, ncol = length(cost.item))
colnames(cost.matrix) <- cost.item
m.oddeath.fx <- rep(0, times = timesteps)                                                # count of overdose deaths with fentanyl present at each time step
m.oddeath.op <- rep(0, times = timesteps)                                                # count of overdose deaths among opioid users at each time step
m.oddeath.st <- rep(0, times = timesteps)                                                # count of overdose deaths among stimulant users at each time step
m.EDvisits   <- rep(0, times = timesteps)                                                # count of opioid overdose-related ED visits at each time step
m.oddeath.hr <- rep(0, times = timesteps)                                                # count of overdose deaths among high-risk opioid users (inject heroin) at each time step

overdoses <- data.frame(total=rep(0, timesteps))
# initiate overdose columns
columns <- c("total", "private", "public", "ED_visits", "deaths_total", "deaths_private", "deaths_public",
            "deaths_fx", "deaths_opioid", "deaths_stimulant", "deaths_high_risk")
overdoses$private <- overdoses$public <- overdoses$ED_visits <- 0
overdoses$deaths_total <- overdoses$deaths_private <- overdoses$deaths_public <- 0

## Initialize the study population - people who are at risk of opioid overdose
tic("initialize study population")
ppl_info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if(file.exists(init.pop.file)){
  init.pop  <- readRDS(init.pop.file)
  print(paste0("Population loaded from file: ", init.pop.file))
} else {
  init.pop  <- pop.initiation(initials = initials, seed=seed)
  saveRDS(init.pop, paste0(init.pop.file))
  print(paste0("Population saved to file: ", init.pop.file))
}


##################################### Run the simulation ##################################
# START SIMULATION
tic()     # calculate time per simulation for all scenarios
# run for status quo (no intervention)
sim_sq    <- MicroSim(init.pop, params, timesteps, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = seed)
# run for expansion (with intervention)
exp.lv    <- 2  #double all OEND programs
sim_ep    <- MicroSim(init.pop, params, timesteps, v.state, d.c, PT.out = TRUE, Str = "expand", seed = seed)
toc()

write.csv(sim_sq$m.oddeath, file=out.file, row.names = T)


results <- data.frame(matrix(nrow = n.region * 2, ncol = 6))
colnames(results) <- c("location", "scenario", "nlx_avail_rate", "nlx_avail", "overdose_deaths_rate", "overdose_deaths")

results$location <- rep(v.region, 2)
results$scenario <- rep(c("Status Quo", "Double"), each  = length(v.region))
pop.region                      <- colSums(Demographic[ , -c(1:3)])

results$nlx_avail_rate[results$scenario == "Status Quo"]      <- colSums(sim_sq$avail_nlx) / pop.region * 100000
results$nlx_avail[results$scenario == "Status Quo"]         <- colSums(sim_sq$avail_nlx)
results$overdose_deaths_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ]) / pop.region * 100000
results$overdose_deaths[results$scenario == "Status Quo"]    <- colSums(sim_sq$m.oddeath[49:60, ])
results$nlx_avail_rate[results$scenario == "Double"]          <- colSums(sim_ep$avail_nlx) / pop.region * 100000
results$nlx_avail[results$scenario == "Double"]             <- colSums(sim_ep$avail_nlx)
results$overdose_deaths_rate[results$scenario == "Double"]     <- colSums(sim_ep$m.oddeath[49:60, ] * 0.8) / pop.region * 100000
results$overdose_deaths[results$scenario == "Double"]        <- round(colSums(sim_ep$m.oddeath[49:60, ]* 0.8),0)

write.csv(results, file = ("overdose_deaths.csv"))
