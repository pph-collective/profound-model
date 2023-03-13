#!/usr/bin/env Rscript

###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2022 #########################
###############################################################################################
# Main module for the microsimulation of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH; Sam Bessy, MSc
# Marshall Lab, Department of Epidemiology, Brown University
# Created: May 06, 2020
# Setting: Massachusetts
###############################################################################################

###############################################################################################
####    Individual-based microsimulation                                                   ####
####    7 health states: prescribed (preb), illicit-injection (il.hr) - high-risk          ####
####                     illicit-noninjection (il.lr) - low-risk                           ####
####                     inactive (inact),  non-opioid drug use (NODU) - stimulant,        ####
####                     relapsed (relap), death (dead)                                    ####
####    1 health event:  Overdose                                                          ####
####    Attributes:      age, sex, residence, race,                                        ####
####                     current state (curr.state), opioid use state (OU.state),          ####
####                     initial state (init.state), initial age (inits.age),              ####
####                     fenatneyl exposure (fx), ever overdosed (ever.od)                 ####
####    Built to inform Naloxone distribution strategies to prevent overdose death         ####
###############################################################################################

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list = ls())
library(argparser)
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(tictoc)
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")

yr_start <- 2016 # starting year of simulation
yr_end <- 2019 # end year of simulation (also the year for evaluation)
discount.rate <- 0.03 # discounting of costs by 3%
# annual.entry <- 16278
# n.add.pop <- round(annual.entry/12, 0) # number of individuals entering the model every month

args <- arg_parser("arguments")
args <- add_argument(args, "--seed", help = "seed for random numbers", default = 2022)
args <- add_argument(args, "--regional", help = "flag to run regional model", flag = TRUE)
args <- add_argument(args, "--outfile", help = "file to store outputs", default = "OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--ppl", help = "file with initial ppl info", default = "Inputs/init_pop.rds")
argv <- parse_args(args)
seed <- as.integer(argv$seed)
init_ppl.file <- argv$ppl
initials$ini.inactive['other'] <- initials$ini.inactive['white'] #TEMPORARY
source("io_setup.R")

Target <- read.xlsx(WB, sheet = "Target")
calib.comp <- Target
calib.comp <- subset(calib.comp, select = -other)
calib.comp <- subset(calib.comp, year < 2020)
calib.comp$white_sim <- rep(0, 12)
calib.comp$black_sim <- rep(0, 12)
calib.comp$hisp_sim  <- rep(0, 12)
calib.comp$wbh_sim   <- rep(0, 12)

## Run simulation ##
for (ss in 3:3){
  seed <- sim.seed[ss]
  params <- sim.data.ls[[ss]]
  params$gw.fx["white"] <- 0.0005
  # params$ini.oud.fx["white"] <- 0.6
  params$nlx.adj <- 1.2
  sim_sq <- MicroSim(init_ppl, params = params, timesteps, agent_states, discount.rate, PT.out = F, strategy = "SQ", seed = seed)        # run for status quo
  calib.comp$white_sim[calib.comp$par == "ODdeaths"] <- colSums(matrix(rowSums(sim_sq$m.oddeath.white), nrow = 12))[1:4]
  calib.comp$black_sim[calib.comp$par == "ODdeaths"] <- colSums(matrix(rowSums(sim_sq$m.oddeath.black), nrow = 12))[1:4]
  calib.comp$hisp_sim[calib.comp$par == "ODdeaths"]  <- colSums(matrix(rowSums(sim_sq$m.oddeath.hisp), nrow = 12))[1:4]
  calib.comp$white_sim[calib.comp$par == "Fx_ODD"] <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "white"], nrow = 12))[1:4] / calib.comp$white_sim[calib.comp$par == "ODdeaths"]
  calib.comp$black_sim[calib.comp$par == "Fx_ODD"] <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "black"], nrow = 12))[1:4] / calib.comp$black_sim[calib.comp$par == "ODdeaths"]
  calib.comp$hisp_sim[calib.comp$par == "Fx_ODD"]  <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "hisp"], nrow = 12))[1:4] / calib.comp$hisp_sim[calib.comp$par == "ODdeaths"]
  calib.comp$white_sim[calib.comp$par == "EDvisits"] <- colSums(matrix(sim_sq$m.EDvisits.race[, "white"], nrow = 12))[1:4]
  calib.comp$black_sim[calib.comp$par == "EDvisits"] <- colSums(matrix(sim_sq$m.EDvisits.race[, "black"], nrow = 12))[1:4]
  calib.comp$hisp_sim[calib.comp$par == "EDvisits"]  <- colSums(matrix(sim_sq$m.EDvisits.race[, "hisp"], nrow = 12))[1:4]
}
calib.comp
