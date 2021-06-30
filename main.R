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
args = commandArgs(trailingOnly=TRUE)
init.pop.file <- "Inputs/InitialPopulation.rds"
## Model setup parameters ##
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"
if (length(args) > 0){
  sw.EMS.ODloc <- args[1]
  seed <- strtoi(args[2])
}

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


# INPUT PARAMETERS
t_start    <- 2016
t_end     <- 2020
ppl_info    <- c("sex", "race", "age", "residence", "curr.state",
                 "OU.state", "init.age", "init.state", "ever.od", "fx")            # information for each model individual
agent_states     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
num_states     <- length(agent_states)                                                     # number of states
num_years        <- t_end-t_start+1
timesteps         <- 12 * num_years                                                           # number of time cycles (in month)
n.rgn       <- length(v.rgn)                                                       # number of regions

# OUTPUT matrices and vectors
# TO_REVIEW: why are these separate vectors instead of one output data frame? And why both for od_death
v.od        <- rep(0, times = timesteps)                                                 # count of overdose events at each time step
v.oddeath   <- rep(0, times = timesteps)                                                 # count of overdose deaths at each time step
m.oddeath   <- matrix(0, nrow = timesteps, ncol = n.rgn)
colnames(m.oddeath) <- v.rgn
v.odpriv    <- rep(0, times = timesteps)                                                 # count of overdose events occurred at private setting at each time step
v.odpubl    <- rep(0, times = timesteps)                                                 # count of overdose events occurred at public setting at each time step
v.deathpriv <- rep(0, times = timesteps)                                                 # count of overdose deaths occurred at private setting at each time step
v.deathpubl <- rep(0, times = timesteps)                                                 # count of overdose deaths occurred at public setting at each time step
v.str       <- c("SQ", "Expand100")                                                # store the strategy names
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
columns <- c("total", "private", "public", "ED_visits", "deaths_total", "deaths_private", "deaths_public"
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
  saveRDS(init.pop, paste0("Inputs/InitialPopulation.rds"))
}
toc()


##################################### Run the simulation ##################################
# START SIMULATION
tic("Simulations: ")
sim_sq    <- MicroSim(init.pop, vparameters, timesteps, agent_states, d.c, PT.out = TRUE, Str = "SQ", seed = seed)  # status quo
sim_ep    <- MicroSim(init.pop, vparameters, timesteps, agent_states, d.c, PT.out = TRUE, Str = "Expand", seed = seed)  # expanded treatment
toc()

write.csv(sim_sq$m.oddeath, file=out.file, row.names = T)


preliminary.results <- data.frame(matrix(nrow = n.rgn * 2, ncol = 6))
colnames(preliminary.results) <- c("location", "scenario", "Rate_Nx", "N_Nx", "Rate_ODdeath", "N_ODdeath")

# TO_REVIEW: What is v.rgn? Why are results from the main script preliminary?
preliminary.results$location <- rep(v.rgn, 2)
preliminary.results$scenario <- rep(c("Status Quo", "Double"), each  = length(v.rgn))
pop.rgn                      <- colSums(Demographic[ , -c(1:3)])
preliminary.results$Rate_Nx[preliminary.results$scenario == "Status Quo"]      <- colSums(sim_sq$n.nlx.mx.str) / pop.rgn * 100000
preliminary.results$N_Nx[preliminary.results$scenario == "Status Quo"]         <- colSums(sim_sq$n.nlx.mx.str)
preliminary.results$Rate_ODdeath[preliminary.results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ]) / pop.rgn * 100000
preliminary.results$N_ODdeath[preliminary.results$scenario == "Status Quo"]    <- colSums(sim_sq$m.oddeath[49:60, ])
preliminary.results$Rate_Nx[preliminary.results$scenario == "Double"]      <- colSums(sim_ep$n.nlx.mx.str) / pop.rgn * 100000
preliminary.results$N_Nx[preliminary.results$scenario == "Double"]         <- colSums(sim_ep$n.nlx.mx.str)
preliminary.results$Rate_ODdeath[preliminary.results$scenario == "Double"] <- colSums(sim_ep$m.oddeath[49:60, ] * 0.8) / pop.rgn * 100000
preliminary.results$N_ODdeath[preliminary.results$scenario == "Double"]    <- round(colSums(sim_ep$m.oddeath[49:60, ]* 0.8),0)

write.csv(preliminary.results, file = ("preliminary_results.csv"))

sum(sim_sq$v.oddeath)
sim_ep$v.oddeath
