#!/usr/bin/env Rscript

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
rm(list = ls())
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

yr_start <- 2016 # starting year of simulation
yr_end <- 2020 # end year of simulation (also the year for evaluation)
d.c <- 0.03 # discounting of costs by 3%
source("io_setup.R")

args <- arg_parser("arguments")
args <- add_argument(args, "--seed", help = "seed for random numbers", default = 2021)
args <- add_argument(args, "--regional", help = "flag to run regional model", flag = TRUE)
args <- add_argument(args, "--outfile", help = "file to store outputs", default = "OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--ppl", help = "file with initial ppl info", default = "Inputs/InitialPopulation.rds")
argv <- parse_args(args)
seed <- as.integer(argv$seed)

init_ppl.file <- argv$ppl
## Model setup parameters ##
sw.EMS.ODloc <- "ov"
out.file <- argv$outfile
if (isTRUE(argv$regional)) {
  sw.EMS.ODloc <- "sp"
}

## Initialize the study population - people who are at risk of opioid overdose
ppl_info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if (file.exists(init_ppl.file)) {
  init_ppl <- readRDS(init_ppl.file)
  print(paste0("Population loaded from file: ", init_ppl.file))
} else {
  init_ppl <- initiate_ppl(initials = initials, seed = seed)
  saveRDS(init_ppl, init_ppl.file)
  print(paste0("Population saved to file: ", init_ppl.file))
}


##################################### Run the simulation ##################################
# START SIMULATION
tic("Simulation time: ") # calculate simulation time
# run for status quo (no intervention)
sim_sq <- MicroSim(init_ppl, params, timesteps, agent_states, d.c, PT.out = TRUE, strategy = "SQ", seed = seed)
# run for expansion (with intervention)
exp.lv <- 2 # double all OEND programs
sim_ep <- MicroSim(init_ppl, params, timesteps, agent_states, d.c, PT.out = TRUE, strategy = "expand", seed = seed)
toc()

write.csv(sim_sq$m.oddeath, file = out.file, row.names = T)


results <- data.frame(matrix(nrow = num_regions * 2, ncol = 6))
colnames(results) <- c("location", "scenario", "nlx_avail_rate", "nlx_avail", "overdose_deaths_rate", "overdose_deaths")

results$location <- rep(v.region, 2)
results$scenario <- rep(c("Status Quo", "Double"), each = length(v.region))
ppl_region <- colSums(Demographic[, -c(1:3)])

results$nlx_avail_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$avail_nlx) / ppl_region * 100000
results$nlx_avail[results$scenario == "Status Quo"] <- colSums(sim_sq$avail_nlx)
results$overdose_deaths_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ]) / ppl_region * 100000
results$overdose_deaths[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ])
results$nlx_avail_rate[results$scenario == "Double"] <- colSums(sim_ep$avail_nlx) / ppl_region * 100000
results$nlx_avail[results$scenario == "Double"] <- colSums(sim_ep$avail_nlx)
results$overdose_deaths_rate[results$scenario == "Double"] <- colSums(sim_ep$m.oddeath[49:60, ] * 0.8) / ppl_region * 100000
results$overdose_deaths[results$scenario == "Double"] <- round(colSums(sim_ep$m.oddeath[49:60, ] * 0.8), 0)

write.csv(results, file = ("overdose_deaths.csv"))
