#!/usr/bin/env Rscript

###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2022 #########################
###############################################################################################
# Main module for the microsimulation of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH; Sam Bessy, MSc
# Marshall Lab, Department of Epidemiology, Brown University
# Created: May 06, 2020
# ATTN: currently not used but will incorporate other main analysis modules here later
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
yr_end <- 2022 # end year of simulation (also the year for evaluation)
discount.rate <- 0.03 # discounting of costs by 3%

args <- arg_parser("arguments")
args <- add_argument(args, "--seed", help = "seed for random numbers", default = 2021)
args <- add_argument(args, "--regional", help = "flag to run regional model", flag = TRUE)
args <- add_argument(args, "--outfile", help = "file to store outputs", default = "OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--ppl", help = "file with initial ppl info", default = "Inputs/init_pop.rds")
argv <- parse_args(args)
seed <- as.integer(argv$seed)
init_ppl.file <- argv$ppl
source("io_setup.R")

##################################### Run simulation ######################################################
############## RI modeling analysis: distributing 10,000 kits through different programs ##################
###########################################################################################################
# Define different program scenarios for distributing additional 10,000 kits and initialize matrices and arrrays for results
sim.seed <- sim.seed[1:500]
scenario.name <- c("Status Quo", "SSP_10K", "MailEvent_10K", "Healthcare_10K")
mx.od.death.last  <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
mx.costs.totl <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
colnames(mx.od.death.last) <- colnames(mx.costs.totl) <-scenario.name
mx.od.death.1st <- mx.od.death.2nd <- mx.od.death.totl <- mx.od.death.wtns.1st <- mx.od.death.wtns.2nd <- mx.od.death.wtns.last <- mx.od.death.wtns.totl <- mx.od.death.last 

array.od.death.rgn <- array(0, dim = c(num_regions, length(scenario.name), length(sim.seed)))
array.oend.nlx.rgn <- array(0, dim = c(num_regions, length(scenario.name), length(sim.seed)))
dimnames(array.od.death.rgn)[[1]] <- dimnames(array.oend.nlx.rgn)[[1]] <- v.region
dimnames(array.od.death.rgn)[[2]] <- dimnames(array.oend.nlx.rgn)[[2]] <- scenario.name
# Simulation based on multiple parameter sets and seeds (calibrated)
for (ss in 1:length(sim.seed)){
  # tic()
  print(paste0("Parameter set: ", ss))
  vparameters.temp <- sim.data.ls[[ss]]
  sim_sq <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = sim.seed[ss])        # run for status quo
  mx.od.death.1st[ss, "Status Quo"]  <- sum(sim_sq$m.oddeath[(timesteps-35):(timesteps-24), ])
  mx.od.death.2nd[ss, "Status Quo"]  <- sum(sim_sq$m.oddeath[(timesteps-23):(timesteps-12), ])
  mx.od.death.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  mx.od.death.totl[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-35):timesteps, ])
  mx.od.death.wtns.1st[ss, "Status Quo"]  <- sum(sim_sq$v.oddeath.w[(timesteps-35):(timesteps-24)])
  mx.od.death.wtns.2nd[ss, "Status Quo"]  <- sum(sim_sq$v.oddeath.w[(timesteps-23):(timesteps-12)])
  mx.od.death.wtns.last[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-11):timesteps])
  mx.od.death.wtns.totl[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-35):timesteps])
  mx.costs.totl[ss, "Status Quo"] <- sim_sq$total.cost
  array.od.death.rgn[ , "Status Quo", ss] <- colSums(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  array.oend.nlx.rgn[ , "Status Quo", ss] <- as.vector(t(sim_sq$n.nlx.OEND.str))
  
  for (jj in 2:length(scenario.name)){
    vparameters.temp$expand.kits <- 10000
    sim_pg <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = scenario.name[jj], seed = sim.seed[ss]) # run for program scenario
    mx.od.death.1st[ss, jj]  <- sum(sim_pg$m.oddeath[(timesteps-35):(timesteps-24), ])
    mx.od.death.2nd[ss, jj]  <- sum(sim_pg$m.oddeath[(timesteps-23):(timesteps-12), ])
    mx.od.death.last[ss, jj] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    mx.od.death.totl[ss, jj] <- sum(sim_pg$m.oddeath[(timesteps-35):timesteps, ])
    mx.od.death.wtns.1st[ss, jj]  <- sum(sim_pg$v.oddeath.w[(timesteps-35):(timesteps-24)])
    mx.od.death.wtns.2nd[ss, jj]  <- sum(sim_pg$v.oddeath.w[(timesteps-23):(timesteps-12)])
    mx.od.death.wtns.last[ss, jj] <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
    mx.od.death.wtns.totl[ss, jj] <- sum(sim_pg$v.oddeath.w[(timesteps-35):timesteps])
    mx.costs.totl[ss, jj] <- sim_pg$total.cost
    array.od.death.rgn[ , jj, ss] <- colSums(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    array.oend.nlx.rgn[ , jj, ss] <- as.vector(t(sim_pg$n.nlx.OEND.str))
  }
  # toc()
}

rgn.results <- data.frame(location = rep(v.region, 4), scenario = rep(scenario.name, each=num_regions), 
                          N_Nx_mean = numeric(num_regions * 4), N_Nx_lower = numeric(num_regions * 4), N_Nx_upper = numeric(num_regions * 4),
                          N_ODDeath_mean = numeric(num_regions * 4), N_ODDeath_lower = numeric(num_regions * 4), N_ODDeath_upper = numeric(num_regions * 4),
                          population = rep(colSums(Demographic[,-c(1:3)]), 4))

for (ii in 1:length(scenario.name)){
  rgn.results$N_Nx_mean[rgn.results$scenario== scenario.name[ii]]  <- apply(array.oend.nlx.rgn, c(1,2), mean)[,scenario.name[ii]]
  rgn.results$N_Nx_lower[rgn.results$scenario== scenario.name[ii]] <- apply(array.oend.nlx.rgn, c(1,2), function(x) quantile(x, probs = c(0.025)))[,scenario.name[ii]]
  rgn.results$N_Nx_upper[rgn.results$scenario== scenario.name[ii]] <- apply(array.oend.nlx.rgn, c(1,2), function(x) quantile(x, probs = c(0.975)))[,scenario.name[ii]]
  rgn.results$N_ODDeath_mean[rgn.results$scenario== scenario.name[ii]]  <- apply(array.od.death.rgn, c(1,2), mean)[,scenario.name[ii]]
  rgn.results$N_ODDeath_lower[rgn.results$scenario== scenario.name[ii]] <- apply(array.od.death.rgn, c(1,2), function(x) quantile(x, probs = c(0.025)))[,scenario.name[ii]]
  rgn.results$N_ODDeath_upper[rgn.results$scenario== scenario.name[ii]] <- apply(array.od.death.rgn, c(1,2), function(x) quantile(x, probs = c(0.975)))[,scenario.name[ii]]
}

detach("package:openxlsx", unload = TRUE)
library(xlsx)

write.xlsx(mx.od.death.1st, file = ("Outputs/RI10K/ODdeaths.xlsx"), sheetName = "1st year", row.names = F)
write.xlsx(mx.od.death.2nd, file = ("Outputs/RI10K/ODdeaths.xlsx"), sheetName = "2nd year", append = T, row.names = F)
write.xlsx(mx.od.death.last, file = ("Outputs/RI10K/ODdeaths.xlsx"), sheetName= "last year", append = T, row.names = F)
write.xlsx(mx.od.death.totl, file = ("Outputs/RI10K/ODdeaths.xlsx"), sheetName= "Total", append = T, row.names = F)
write.xlsx(mx.od.death.wtns.1st, file = ("Outputs/RI10K/WitnessedDeaths.xlsx"), sheetName = "1st year", row.names = F)
write.xlsx(mx.od.death.wtns.2nd, file = ("Outputs/RI10K/WitnessedDeaths.xlsx"), sheetName = "2nd year", append = T, row.names = F)
write.xlsx(mx.od.death.wtns.last, file = ("Outputs/RI10K/WitnessedDeaths.xlsx"), sheetName= "last year", append = T, row.names = F)
write.xlsx(mx.od.death.wtns.totl, file = ("Outputs/RI10K/WitnessedDeaths.xlsx"), sheetName= "Total", append = T, row.names = F)
write.xlsx(mx.costs.totl, file = ("Outputs/RI10K/Total_costs.xlsx"), row.names = F)
write.csv(rgn.results, file = ("Outputs/RI10K/Results_bytown.csv"), row.names = F)
