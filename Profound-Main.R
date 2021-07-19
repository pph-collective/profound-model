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
rm(list=ls())
args = commandArgs(trailingOnly=TRUE)

## Model setup parameters ##
seed         <- 2021
sw.EMS.ODloc <- "ov"  #Please choose from "ov" (using average overall) or "sp" (region-specific) for overdose setting parameter, default is "ov"

# install.packages("rstudioapi")
# library(rstudioapi)
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
source("Profound-Function-PopInitialization.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-Function-Microsimulation.R")
source("Profound-DecisionTree.R")
source("Profound-DataInput.R")
source("Profound-Function-NxAvailAlgm.R")
source("Profound-CEA.R")


# INPUT PARAMETERS
yr.first    <- 2016        # starting year of simulation 
yr.last     <- 2020        # end year of simulation (also the year for evaluation)
d.c         <- 0.03        # discounting of costs by 3%                                    

source("Profound-InputOutput-Setup.R")

## Initialize the study population - people who are at risk of opioid overdose
pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if(file.exists("Inputs/InitialPopulation.rds")){
  init.pop  <- readRDS("Inputs/InitialPopulation.rds")
  print(paste0("Population loaded from file: ", "Inputs/InitialPopulation.rds"))
} else {
  tic()
  init.pop  <- pop.initiation(initials = initials, seed=seed)
  saveRDS(init.pop, paste0("Inputs/InitialPopulation.rds"))
  toc()
}


##################################### Run the simulation ##################################
# START SIMULATION
tic()     # calculate time per simulation for all scenarios
# run for status quo (no intervention)
sim_sq    <- MicroSim(init.pop, vparameters, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = seed)
# run for expansion (with intervention)
exp.lv    <- 2  #double all OEND programs
sim_ep    <- MicroSim(init.pop, vparameters, n.t, v.state, d.c, PT.out = TRUE, Str = "expand", seed = seed)
toc()