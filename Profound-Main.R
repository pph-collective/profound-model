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
args <- commandArgs(trailingOnly = TRUE)

## Model setup parameters ##
seed <- 2021
sw.EMS.ODloc <- "ov" # Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

# install.packages("rstudioapi")
# library(rstudioapi)
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
source("population.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-Function-Microsimulation.R")
source("Profound-DecisionTree.R")
source("Profound-DataInput.R")
source("Profound-Function-NxAvailAlgm.R")
source("Profound-CEA.R")


# INPUT PARAMETERS
yr.first <- 2016 # starting year of simulation
yr.last <- 2020 # end year of simulation (also the year for evaluation)
d.c <- 0.03 # discounting of costs by 3%

source("io_setup.R")

## Initialize the study population - people who are at risk of opioid overdose
pop.info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if (file.exists("Inputs/InitialPopulation.rds")) {
  init_ppl <- readRDS("Inputs/InitialPopulation.rds")
  print(paste0("Population loaded from file: ", "Inputs/InitialPopulation.rds"))
} else {
  tic()
  init_ppl <- pop.initiation(initials = initials, seed = seed)
  saveRDS(init_ppl, paste0("Inputs/InitialPopulation.rds"))
  toc()
}


##################################### Run the simulation ##################################
# START SIMULATION
tic() # calculate time per simulation for all scenarios
# run for status quo (no intervention)
sim_sq <- MicroSim(init_ppl, params, timesteps, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = seed)
# run for expansion (with intervention)
exp.lv <- 2 # double all OEND programs
sim_ep <- MicroSim(init_ppl, params, timesteps, v.state, d.c, PT.out = TRUE, Str = "expand", seed = seed)
toc()

write.csv(sim_sq$m.oddeath, file = "OverdoseDeath_RIV1.0.csv", row.names = T)


preliminary.results <- data.frame(matrix(nrow = n.rgn * 2, ncol = 6))
x <- c("location", "scenario", "Rate_Nx", "N_Nx", "Rate_ODdeath", "N_ODdeath")
colnames(preliminary.results) <- x

preliminary.results$location <- rep(v.rgn, 2)
preliminary.results$scenario <- rep(c("Status Quo", "Double"), each = length(v.rgn))
pop.rgn <- colSums(Demographic[, -c(1:3)])
preliminary.results$Rate_Nx[preliminary.results$scenario == "Status Quo"] <- colSums(sim_sq$n.nlx.mx.str) / pop.rgn * 100000
preliminary.results$N_Nx[preliminary.results$scenario == "Status Quo"] <- colSums(sim_sq$n.nlx.mx.str)
preliminary.results$Rate_ODdeath[preliminary.results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ]) / pop.rgn * 100000
preliminary.results$N_ODdeath[preliminary.results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ])
preliminary.results$Rate_Nx[preliminary.results$scenario == "Double"] <- colSums(sim_ep$n.nlx.mx.str) / pop.rgn * 100000
preliminary.results$N_Nx[preliminary.results$scenario == "Double"] <- colSums(sim_ep$n.nlx.mx.str)
preliminary.results$Rate_ODdeath[preliminary.results$scenario == "Double"] <- colSums(sim_ep$m.oddeath[49:60, ] * 0.8) / pop.rgn * 100000
preliminary.results$N_ODdeath[preliminary.results$scenario == "Double"] <- round(colSums(sim_ep$m.oddeath[49:60, ] * 0.8), 0)

write.csv(preliminary.results, file = ("preliminary.results.csv"))
