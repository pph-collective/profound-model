###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2020 #########################
###############################################################################################
# Main module for the microsimulation of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: May 06, 2020
# Last update: May 28, 2020
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

library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(here)

# parse command line args
args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
  WB.path <- here("..", "Inputs", "MasterTable.xlsx")
} else {
  WB.path <- args[1]
}

source(here("Profound-Function-PopInitialization.R"))
source(here("Profound-Function-TransitionProbability.R"))
source(here("Profound-Function-Microsimulation.R"))
source(here("Profound-DecisionTree.R"))
source(here("Profound-DataInput.R"))
source(here("Profound-Function-NxAvailAlgm.R"))
source(here("Profound-CEA.R"))


# INPUT PARAMETERS
yr.first    <- 2016
yr.last     <- 2020
pop.info    <- c("sex", "race", "age", "residence", "curr.state",
                 "OU.state", "init.age", "init.state", "ever.od", "fx")            # information for each model individual
v.state     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
n.state     <- length(v.state)                                                     # number of states
n.yr        <- yr.last-yr.first+1
n.t         <- 12 * n.yr                                                           # number of time cycles (in month)
n.rgn       <- length(v.rgn)                                                       # number of regions
v.od        <- rep(0, times = n.t)                                                 # count of overdose events at each time step
v.oddeath   <- rep(0, times = n.t)                                                 # count of overdose deaths at each time step
m.oddeath   <- matrix(0, nrow = n.t, ncol = n.rgn)
colnames(m.oddeath) <- v.rgn
v.str       <- c("SQ", "Expand50", "All+200", "R2+600")                            # store the strategy names
d.c         <- 0.03                                                                # discounting of costs by 3%
cost.item   <- c("TotalCost", "NxCost")
cost.matrix <- matrix(0, nrow=n.t, ncol = length(cost.item))
colnames(cost.matrix) <- cost.item
array.Nx    <- array.Nx.full[dimnames(array.Nx.full)[[1]]>=yr.first, , ]
init.Nx     <- array.Nx.full[dimnames(array.Nx.full)[[1]]==yr.first-1, , ]


# # Initialize the study population - people who are at risk of opioid overdose
pop.info  <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
init.pop.file = here("..", "Inputs", "InitialPopulation.rds")
if(file.exists(init.pop.file)){
  init.pop  <- readRDS(init.pop.file)
} else {
  tic()
  init.pop  <- pop.initiation(initials = initials, seed=2021)
  saveRDS(init.pop, init.pop.file)
  toc()
}

# Overdose probability matrix (per month)
od.matrix             <- matrix(0, nrow = 4, ncol = 2)
rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
colnames(od.matrix)   <- c("first", "subs")
od.matrix["preb", "subs"]   <- od.preb.sub
od.matrix["il.lr", "subs"]  <- od.il.lr.sub
od.matrix["il.hr", "subs"]  <- od.il.lr.sub * multi.hr
od.matrix["NODU", "subs"]   <- od.NODU.sub
od.matrix[ , "first"]       <- od.matrix[ , "subs"] / multi.sub

# Baseline mortality excluding overdose (per month)
mor.matrix                  <- matrix(0, nrow = 2, ncol = length(mor.gp))
rownames(mor.matrix)        <- c("bg", "drug")
colnames(mor.matrix)        <- mor.gp
mor.matrix["bg", ]          <- mor.bg
mor.matrix["drug", ]        <- mor.drug


# # Naloxone distribution
# n.nlx.v       <- c(0, 1000, 200)
# n.od_death.v  <- c(150, 800, 250)
# nlx.adj       <- 1                  # adjuster for naloxone availability given the ratio between naloxone kits to od_deaths in a region
#
# # Spatial dynamics to access naloxone kits
# R1.nlx  <-  c(0.7, 0.2, 0.1)
# R2.nlx  <-  c(0.15, 0.8, 0.05)
# R3.nlx  <-  c(0.1, 0.25, 0.65)
# acc.nlx.matrix <- rbind(R1.nlx, R2.nlx, R3.nlx)
# rownames(acc.nlx.matrix)  <- c("R1.from", "R2.from", "R3.from")
# colnames(acc.nlx.matrix)  <- c("R1.to", "R2.to", "R3.to")


##################################### Run the simulation ##################################
# START SIMULATION
p = Sys.time()
sim_sq   <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = 100)        # run for no treatment
# sim_ep50 <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "Expand50", seed = 100)  # run for treatment
# sim_a200 <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "All+200", seed = 100)   # run for treatment
# sim_r2   <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "R2+600", seed = 100)    # run for treatment

comp.time = Sys.time() - p

comp.time

write.csv(sim_sq$m.oddeath, file="OverdoseDeath_RIV1.0.csv", row.names = T)

sum(sim_sq$v.oddeath)
# sum(sim_ep50$v.oddeath)
# sum(sim_a200$v.oddeath)
# sum(sim_r2$v.oddeath)
