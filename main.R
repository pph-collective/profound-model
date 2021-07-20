#!/usr/bin/env Rscript

###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2020 #########################
###############################################################################################
# Main module for the microsimulation of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH; Sam Bessy, MSc
# Marshall Lab, Department of Epidemiology, Brown University
# Created: May 06, 2020
# Last update: July 17, 2021
# ATTN: currently not used but will incorporate other main analysis modules here later
###############################################################################################

###############################################################################################
####    Individual-based microsimulation                                                   ####
####    7 health states: prescribed (preb), unregulated-injection (unreg.inj)              ####
####                     unregulated-noninjection (unreg.nin)                              ####
####                     inactive (inact),  non-opioid drug use (NODU) - stimulant,        ####
####                     relapsed (relap), death (dead)                                    ####
####    1 health event:  Overdose                                                          ####
####    Attributes:      age, sex, residence, race,                                        ####
####                     current state (curr.state), opioid use state (OU.state),          ####
####                     initial state (init.state), initial age (inits.age),              ####
####                     fenatneyl exposure (fx), ever overdosed (ever.od)                 ####
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
sw.EMS.ODloc <- "overall"
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
sim_ep <- MicroSim(init_ppl, params, timesteps, v.state, d.c, PT.out = TRUE, strategy = "expand", seed = seed)
toc()
