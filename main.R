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
####    Built to inform Naloxone distribution strategies to prevent overdose death         ####
###############################################################################################

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list = ls())
library("argparser")
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
scenario.name <- c("Status Quo", "SSP_10K", "StreetOutreach_10K", "MailEvent_10K", "Healthcare_10K")
mx.od.death.last <- mx.od.death.totl <- mx.od.death.wtns.last <- mx.od.death.wtns.totl <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
mx.costs.totl <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
colnames(mx.od.death.last) <- colnames(mx.od.death.totl) <- colnames(mx.od.death.wtns.last) <- colnames(mx.od.death.wtns.totl) <- colnames(mx.costs.totl) <-scenario.name
array.od.death.rgn <- array(0, dim = c(num_regions, length(scenario.name), length(sim.seed)))
array.oend.nlx.rgn <- array(0, dim = c(num_regions, length(scenario.name), length(sim.seed)))
dimnames(array.od.death.rgn)[[1]] <- dimnames(array.oend.nlx.rgn)[[1]] <- v.region
dimnames(array.od.death.rgn)[[2]] <- dimnames(array.oend.nlx.rgn)[[2]] <- scenario.name
# Simulation based on multiple parameter sets and seeds (calibrated)
for (ss in 1:length(sim.seed)){
  # tic()
  print(paste0("Parameter set: ", ss))
  vparameters.temp <- sim.data.ls[[ss]]
  vparameters.temp$NxDataPharm$pe <- round(vparameters.temp$NxDataPharm$pe / 2.69, 0) # ATTN: This is temporary, used to manually adjust pharmacy naloxone
  # vparameters.temp$gw.m.2inact <- 0.0059  # ATTN: This is temporary, used to manually add monthly growth for transition into inactive (for MAT increase)
  vparameters.temp$gw.m.2inact <- 0  # ATTN: This is temporary, used to manually add monthly growth for transition into inactive (for MAT increase)
  vparameters.temp$nlx.adj <- 1.2  # ATTN: This is temporary, used to manually adjust naloxone adjustment term to 1
  vparameters.temp$mor_bl <- 0.0588 * 0.58
  vparameters.temp$mor_nx <- 0.00899 * 0.58
  vparameters.temp$OD_cess <- 0
  vparameters.temp$gw.fx <- 0.05
  vparameters.temp$multi.relap <- 3
  vparameters.temp$multi.fx <- 6
  # vparameters.temp$p.preb2inact.ini <- vparameters.temp$p.preb2inact
  vparameters.temp$p.preb2inact.ini <- 0.005
  vparameters.temp[(names(vparameters.temp)=='p.preb2inact')] <- NULL # ATTN: This is temporary, used to fix a recent change to 
  vparameters.temp$p.il.lr2inact.ini <- vparameters.temp$p.il.lr2inact
  vparameters.temp[(names(vparameters.temp)=='p.preb2inact')] <- NULL # ATTN: This is temporary, used to manually adjust naloxone adjustment term to 1
  vparameters.temp$p.il.hr2inact.ini <- vparameters.temp$p.il.hr2inact
  vparameters.temp[(names(vparameters.temp)=='p.preb2inact')] <- NULL # ATTN: This is temporary, used to manually adjust naloxone adjustment term to 1
  vparameters.temp$Low2Priv <- 0.3
  vparameters.temp$High2Pub <- 0.1
  # vparameters.temp$p.inact2relap <- 0.004
  # vparameters.temp$overdose_probs = vparameters.temp$overdose_probs * 0.9
    
  # status quo scenario
  sim_sq <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = sim.seed[ss])        # run for status quo
  sim_sq <- MicroSim(init_ppl, params = params, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = 2002)        # run for status quo
  
  mx.od.death.last[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  mx.od.death.totl[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-35):timesteps, ])
  mx.od.death.wtns.last[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-11):timesteps])
  mx.od.death.wtns.totl[ss, "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-35):timesteps])
  mx.costs.totl[ss, "Status Quo"] <- sim_sq$total.cost
  array.od.death.rgn[ , "Status Quo", ss] <- colSums(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  array.oend.nlx.rgn[ , "Status Quo", ss] <- sim_sq$n.nlx.OEND.str
  
  for (jj in 2:length(scenario.name)){
    vparameters.temp$expand.kits <- 10000
    sim_pg <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = scenario.name[jj], seed = sim.seed[ss]) # run for program scenario
    mx.od.death.last[ss, jj] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    mx.od.death.totl[ss, jj] <- sum(sim_pg$m.oddeath[(timesteps-35):timesteps, ])
    mx.od.death.wtns.last[ss, jj] <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
    mx.od.death.wtns.totl[ss, jj] <- sum(sim_pg$v.oddeath.w[(timesteps-35):timesteps])
    mx.costs.totl[ss, jj] <- sim_pg$total.cost
    array.od.death.rgn[ , jj, ss] <- colSums(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    array.oend.nlx.rgn[ , jj, ss] <- sim_pg$n.nlx.OEND.str
  }
  # toc()
}


results.od.death.rgn <- apply(array.od.death.rgn, c(1,2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
results.oend.nlx.rgn <- apply(array.oend.nlx.rgn, c(1,2), function(x) quantile(x, probs = c(0.025, 0.5, 0.975)))
new.column.name <- as.vector(t(outer(dimnames(results.od.death.rgn)[[3]], dimnames(results.od.death.rgn)[[1]], paste, sep="_")))
for(ii in 1:dim(results.od.death.rgn)[3]){
  if (ii == 1){
    mx.results.od.death.rgn <- t(results.od.death.rgn[,,ii])
    mx.results.oend.nlx.rgn <- t(results.oend.nlx.rgn[,,ii])
  } else {
    mx.results.od.death.rgn <- cbind(mx.results.od.death.rgn, t(results.od.death.rgn[,,ii]))
    mx.results.oend.nlx.rgn <- cbind(mx.results.oend.nlx.rgn, t(results.oend.nlx.rgn[,,ii]))
  }
}
colnames(mx.results.od.death.rgn) <- colnames(mx.results.oend.nlx.rgn) <- new.column.name

ppl_region <- colSums(Demographic[, -c(1:3)])
mx.results.od.death.rgn <- cbind(mx.results.od.death.rgn, ppl_region)
mx.results.oend.nlx.rgn <- cbind(mx.results.oend.nlx.rgn, ppl_region)

write.xlsx(mx.od.death.last, file = ("Outputs/RI10K/ODdeaths_LastYear.xlsx"), row.names = F)
write.xlsx(mx.od.death.totl, file = ("Outputs/RI10K/ODdeaths_Total3Y.xlsx"), row.names = F)
write.xlsx(mx.od.death.wtns.last, file = ("Outputs/RI10K/WitnessedDeaths_LastYear.xlsx"), row.names = F)
write.xlsx(mx.od.death.wtns.totl, file = ("Outputs/RI10K/WitnessedDeaths_Total3Y.xlsx"), row.names = F)
write.xlsx(mx.costs.totl, file = ("Outputs/RI10K/Total_costs.xlsx"), row.names = F)
write.xlsx(mx.results.od.death.rgn, file = ("Outputs/RI10K/ODdeaths_bytown.xlsx"), row.names = T)
write.xlsx(mx.results.oend.nlx.rgn, file = ("Outputs/RI10K/OENDnaloxone_bytown.xlsx"), row.names = T)
