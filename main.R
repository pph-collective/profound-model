#!/usr/bin/env Rscript

################################################################################
############### PROFOUND Naloxone Distribution model ### 2020 ##################
################################################################################
# Main module for the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH; Sam Bessy, MSc
# Marshall Lab, Department of Epidemiology, Brown University
# Created: May 06, 2020
# Last update: July 17, 2021
#
################################################################################

################################################################################
####    Individual-based microsimulation                                    ####
####    7 agent states: prescribed (preb),                                  ####
####                     unregulated-injection (unreg.inj)                  ####
####                     unregulated-noninjection (unreg.nin)               ####
####                     inactive (inact),                                  ####
####                     non-opioid drug use (NODU) - stimulant,            ####
####                     relapsed (relap),                                  ####
####                     death (dead)                                       ####
####    1 health event:  Overdose                                           ####
####    Attributes:      age, sex, residence, race,                         ####
####                     current state (curr.state),                        ####
####                     opioid use state (OU.state),                       ####
####                     initial state (init.state),                        ####
####                     initial age (inits.age),                           ####
####                     fenatneyl exposure (fx),                           ####
####                     ever overdosed (ever.od)                           ####
####  Projects effects of naloxone distribution strategies on overdose      ####
################################################################################

################################################################################
# 1. SET directpry and workspace
################################################################################
rm(list = ls())
library("argparser")
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(tictoc)
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")
source("parse_params.R")
source("io_setup.R")
source("data_input.R")

# parse command line arguments--------------------------------------------------
args <- arg_parser("arguments")
# args <- add_argument(args,
#   "--seed",
#   help = "seed for random numbers",
#   default = 2021
# )
# args <- add_argument(args,
#   "--regional",
#   help = "flag to run regional model",
#   flag = TRUE
# )
# args <- add_argument(args,
#   "--outfile",
#   help = "file to store outputs",
#   default = "OverdoseDeath_RIV1_0.csv"
# )
# args <- add_argument(args,
#   "--ppl",
#   help = "file with initial ppl info",
#   default = "Inputs/init_pop.rds"
# )
# args <- add_argument(args,
#   "--cores",
#   help = "number of cores to run in parallel. If none, do not open socket.",
#   default = 1
# )
args <- add_argument(args, "--inputs", help = "file containing input values", default = "params/base_params.yml")
argv <- parse_args(args)

inputs <- parse_inputs(argv$inputs)
d.c <- inputs$discount # discounting of costs by 3%
main_table <- inputs$main_table

sw.EMS.ODloc <- inputs$strat
out.file <- inputs$outfile
# init_ppl.file <- inputs$init_ppl



# add parameters
params <- input_setup(inputs)
# create output table
output <- output_setup(params)
data <- data_input(params)  # empirical
## Initialize the study population - people who are at risk of opioid overdose
ppl_info <- c(
  "sex",
  "race",
  "age",
  "residence",
  "curr.state",
  "OU.state",
  "init.age",
  "init.state",
  "ever.od",
  "fx"
)

if (file.exists(inputs$init_ppl_file)) { # import pop if possible
  init_ppl <- readRDS(inputs$init_ppl_file)
  print(paste0("Population loaded from file: ", inputs$init_ppl_file))
} else { # otherwise, create pop
  init_ppl <- initiate_ppl(data, seed = inputs$seed)
  saveRDS(init_ppl, inputs$init_ppl_file)
  print(paste0("Population saved to file: ", inputs$init_ppl_file))
}


# Run the simulation =============================
# run for status quo (no intervention)
# sim_sq <- MicroSim(init_ppl, params, output, agent_states, d.c, TRUE, "SQ", seed)
# # run for expansion (with intervention)
# exp.lv <- 2 # double all OEND programs
# sim_ep <- MicroSim(init_ppl, params, output, agent_states, d.c, TRUE, "expand", seed)


# results <- data.frame(matrix(nrow = length(v.region) * 2, ncol = 6))
# colnames(results) <- c("location", "scenario", "nlx_avail_rate", "nlx_avail", "overdose_deaths_rate", "overdose_deaths")

# results$location <- rep(v.region, 2)
# results$scenario <- rep(c("Status Quo", "Double"), each = length(v.region))
# ppl_region <- colSums(params$Demographic[, -c(1:3)])

# results$nlx_avail_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$avail_nlx) / ppl_region * 100000
# results$nlx_avail[results$scenario == "Status Quo"] <- colSums(sim_sq$avail_nlx)
# results$overdose_deaths_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ]) / ppl_region * 100000
# results$overdose_deaths[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ])
# results$nlx_avail_rate[results$scenario == "Double"] <- colSums(sim_ep$avail_nlx) / ppl_region * 100000
# results$nlx_avail[results$scenario == "Double"] <- colSums(sim_ep$avail_nlx)
# results$overdose_deaths_rate[results$scenario == "Double"] <- colSums(sim_ep$m.oddeath[49:60, ] * 0.8) / ppl_region * 100000
# results$overdose_deaths[results$scenario == "Double"] <- round(colSums(sim_ep$m.oddeath[49:60, ] * 0.8), 0)
# write.csv(results, file = ("overdose_deaths.csv"))

main <- function(init_ppl, params, data, output, timesteps, agent_states, d.c, PT.out, expansion, seed, scenario){
  # what I want to do: for each scenario, run the simulation. If it's a program, send it to the program eval to change the probs. Otherwise, scale as needed

  sim_sq <- MicroSim(init_ppl, params, data, output, agent_states, d.c, TRUE, "SQ", seed)

}

main(init_ppl, params, data, output, timesteps, agent_states, d.c, PT.out, expansion, seed, scenario)
