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
yr_end <- 2023 # end year of simulation (also the year for evaluation)
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

##################################### Run simulation #####################################################
############## MA modeling analysis: Compare non-targeted VS. equity-focused expansion  ##################
##########################################################################################################
sim.seed <- sim.seed[1]
scenario.name <- c("Status Quo", "NT_40ratio", "NT_60ratio", "NT_80ratio", "EQ_40ratio", "EQ_60ratio", "EQ_80ratio")
ratio.level <-  c(40, 60, 80, 40, 60, 80)
mx.od.death.last  <- matrix(0, nrow = 4, ncol = length(scenario.name))
# mx.costs.totl <- matrix(0, nrow = length(sim.seed), ncol = length(scenario.name))
colnames(mx.od.death.last) <- scenario.name
rownames(mx.od.death.last) <- c("total", "white", "black", "hisp")
mx.od.death.1st <- mx.od.death.2nd <- mx.od.death.totl <- mx.od.death.wtns.1st <- mx.od.death.wtns.2nd <- mx.od.death.wtns.last <- mx.od.death.wtns.totl <- mx.od.death.last 

array.od.death.rgn <- array(0, dim = c(num_regions, 3, length(scenario.name)))
array.oend.nlx.rgn <- array(0, dim = c(num_regions, 3, length(scenario.name)))
dimnames(array.od.death.rgn)[[1]] <- dimnames(array.oend.nlx.rgn)[[1]] <- v.region
dimnames(array.od.death.rgn)[[2]] <- dimnames(array.oend.nlx.rgn)[[2]] <- c("white", "black", "hispanic")
dimnames(array.od.death.rgn)[[3]] <- dimnames(array.oend.nlx.rgn)[[3]] <- scenario.name
vparameters.temp <- sim.data.ls[[3]]
params$overdose_probs <- vparameters.temp$overdose_probs
params$multi.fx <- vparameters.temp$multi.fx
params$multi.relap <- vparameters.temp$multi.relap
params$multi.sub <- vparameters.temp$multi.sub
params$mor_bl <- vparameters.temp$mor_bl
params$mor_nx <- vparameters.temp$mor_nx
params$rr_mor_EMS <- vparameters.temp$rr_mor_EMS
params$OD_wit_pub <- vparameters.temp$OD_wit_pub
params$OD_911_pub <- vparameters.temp$OD_911_pub
params$OD_hosp    <- vparameters.temp$OD_hosp
params$nlx.adj    <- vparameters.temp$nlx.adj*1.2
params$OD_loc_pub <- vparameters.temp$OD_loc_pub
params$OD_wit_priv<- vparameters.temp$OD_wit_priv
params$OD_911_priv<- vparameters.temp$OD_911_priv

Target <- read.xlsx(WB, sheet = "Target")
calib.comp <- Target
calib.comp <- subset(calib.comp, select = -other)
calib.comp$white_sim <- rep(0, 15)
calib.comp$black_sim <- rep(0, 15)
calib.comp$hisp_sim  <- rep(0, 15)
calib.comp$wbh_sim   <- rep(0, 15)

## Hand tune model parameters
params$overdose_probs[1,] <- params$overdose_probs[1,]/2
# params$overdose_probs[2:3,] <- params$overdose_probs[2:3,]/2
# params$overdose_probs[, 1] <- params$overdose_probs[, 2] / 1
params$ini.oud.fx["white"] <- 0.55
params$ini.oud.fx["black"] <- 0.35
params$ini.oud.fx["hisp"]  <- 0.45
params$gw.fx <- c(0, 0.0005, 0.0005)
params$multi.fx <- 8
params$gw.m.2inact <- 0.0059/2
params$p.inact2relap <- 0.0115380920094443
# params$gw.NODU.fx.ab.yr <- 0
params$ini.NOUD.fx["white"] <- 0.1007
params$ini.NOUD.fx["black"] <- 0.1007
params$ini.NOUD.fx["hisp"]  <- 0.1007 
params$covid.NOUD.fx['white'] <- 1
params$covid.NOUD.fx['black'] <-3
params$covid.NOUD.fx['hisp'] <- 1.2
params$covid.rd.2inact.white <- 0
params$covid.rd.2inact.black <- 0.2
params$covid.rd.2inact.hisp  <- 0.32
params$cap <- 0.99
# params$multi.relap <- 1
# params$mortality_probs$mor.drug[] <- params$mortality_probs$mor.drug[] / 2
# params$mortality_probs$mor.drug <- params$mortality_probs$mor.drug

## Hand tune transition probabilities
params$p.preb2il.lr <- params$p.preb2il.lr
params$p.il.lr2il.hr <- params$p.il.lr2il.hr

# params$mortality_probs$mor.drug[] <- 0 
# params$gw.m.2inact  <- 0
# params$OD_cess[] <- 0
# detach("package:openxlsx", unload = TRUE)
# library(xlsx)

seed.sameple <- c(610, 6088, 4488, 803, 2020, 602, 625, 125, 909, 52390)
result.cost <- matrix(0, nrow= length(seed.sameple), ncol = length(scenario.name))
colnames(result.cost) <- scenario.name
result.qaly <- result.cost
result.qaly.w <- result.cost
result.qaly.b <- result.cost
result.qaly.h <- result.cost

## Run simulation ##
for (ss in 1:10){
  seed <- seed.sameple[ss]
  sim_sq <- MicroSim(init_ppl, params = params, timesteps, agent_states, discount.rate, PT.out = F, strategy = "SQ", seed = seed)        # run for status quo

  result.cost[ss, 1] <- sim_sq$total.cost
  result.qaly[ss, 1] <- sim_sq$total.qaly
  result.qaly.w[ss, 1] <- sim_sq$total.qaly.w
  result.qaly.b[ss, 1] <- sim_sq$total.qaly.b
  result.qaly.h[ss, 1] <- sim_sq$total.qaly.h
  # mx.od.death.1st["total", "Status Quo"]  <- sum(sim_sq$m.oddeath[(timesteps-35):(timesteps-24), ])
  # mx.od.death.1st["white", "Status Quo"]  <- sum(sim_sq$m.oddeath.white[(timesteps-35):(timesteps-24), ])
  # mx.od.death.1st["black", "Status Quo"]  <- sum(sim_sq$m.oddeath.black[(timesteps-35):(timesteps-24), ])
  # mx.od.death.1st["hisp", "Status Quo"]   <- sum(sim_sq$m.oddeath.hisp[(timesteps-35):(timesteps-24), ])
  # mx.od.death.2nd["total", "Status Quo"]  <- sum(sim_sq$m.oddeath[(timesteps-23):(timesteps-12), ])
  # mx.od.death.2nd["white", "Status Quo"]  <- sum(sim_sq$m.oddeath.white[(timesteps-23):(timesteps-12), ])
  # mx.od.death.2nd["black", "Status Quo"]  <- sum(sim_sq$m.oddeath.black[(timesteps-23):(timesteps-12), ])
  # mx.od.death.2nd["hisp", "Status Quo"]   <- sum(sim_sq$m.oddeath.hisp[(timesteps-23):(timesteps-12), ])
  # mx.od.death.last["total", "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-11):timesteps, ])
  # mx.od.death.last["white", "Status Quo"] <- sum(sim_sq$m.oddeath.white[(timesteps-11):timesteps, ])
  # mx.od.death.last["black", "Status Quo"] <- sum(sim_sq$m.oddeath.black[(timesteps-11):timesteps, ])
  # mx.od.death.last["hisp", "Status Quo"]  <- sum(sim_sq$m.oddeath.hisp[(timesteps-11):timesteps, ])
  # mx.od.death.totl["total", "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps-35):timesteps, ])
  # mx.od.death.totl["white", "Status Quo"] <- sum(sim_sq$m.oddeath.white[(timesteps-35):timesteps, ])
  # mx.od.death.totl["black", "Status Quo"] <- sum(sim_sq$m.oddeath.black[(timesteps-35):timesteps, ])
  # mx.od.death.totl["hisp", "Status Quo"]  <- sum(sim_sq$m.oddeath.hisp[(timesteps-35):timesteps, ])
  # mx.od.death.wtns.1st["total", "Status Quo"]  <- sum(sim_sq$v.oddeath.w[(timesteps-35):(timesteps-24)])
  # mx.od.death.wtns.1st["white", "Status Quo"]  <- sum(sim_sq$v.oddeath.w.race[(timesteps-35):(timesteps-24), "white"])
  # mx.od.death.wtns.1st["black", "Status Quo"]  <- sum(sim_sq$v.oddeath.w.race[(timesteps-35):(timesteps-24), "black"])
  # mx.od.death.wtns.1st["hisp", "Status Quo"]   <- sum(sim_sq$v.oddeath.w.race[(timesteps-35):(timesteps-24), "hisp"])
  # mx.od.death.wtns.2nd["total", "Status Quo"]  <- sum(sim_sq$v.oddeath.w[(timesteps-23):(timesteps-12)])
  # mx.od.death.wtns.2nd["white", "Status Quo"]  <- sum(sim_sq$v.oddeath.w.race[(timesteps-23):(timesteps-12), "white"])
  # mx.od.death.wtns.2nd["black", "Status Quo"]  <- sum(sim_sq$v.oddeath.w.race[(timesteps-23):(timesteps-12), "black"])
  # mx.od.death.wtns.2nd["hisp", "Status Quo"]   <- sum(sim_sq$v.oddeath.w.race[(timesteps-23):(timesteps-12), "hisp"])
  # mx.od.death.wtns.last["total", "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-11):timesteps])
  # mx.od.death.wtns.last["white", "Status Quo"] <- sum(sim_sq$v.oddeath.w.race[(timesteps-11):timesteps, "white"])
  # mx.od.death.wtns.last["black", "Status Quo"] <- sum(sim_sq$v.oddeath.w.race[(timesteps-11):timesteps, "black"])
  # mx.od.death.wtns.last["hisp", "Status Quo"]  <- sum(sim_sq$v.oddeath.w.race[(timesteps-11):timesteps, "hisp"])
  # mx.od.death.wtns.totl["total", "Status Quo"] <- sum(sim_sq$v.oddeath.w[(timesteps-35):timesteps])
  # mx.od.death.wtns.totl["white", "Status Quo"] <- sum(sim_sq$v.oddeath.w.race[(timesteps-35):timesteps, "white"])
  # mx.od.death.wtns.totl["black", "Status Quo"] <- sum(sim_sq$v.oddeath.w.race[(timesteps-35):timesteps, "black"])
  # mx.od.death.wtns.totl["hisp", "Status Quo"]  <- sum(sim_sq$v.oddeath.w.race[(timesteps-35):timesteps, "hisp"])
  calib.comp$white_sim[calib.comp$par == "ODdeaths"] <- colSums(matrix(rowSums(sim_sq$m.oddeath.white), nrow = 12))[1:5]
  calib.comp$black_sim[calib.comp$par == "ODdeaths"] <- colSums(matrix(rowSums(sim_sq$m.oddeath.black), nrow = 12))[1:5]
  calib.comp$hisp_sim[calib.comp$par == "ODdeaths"]  <- colSums(matrix(rowSums(sim_sq$m.oddeath.hisp), nrow = 12))[1:5]
  calib.comp$white_sim[calib.comp$par == "Fx_ODD"] <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "white"], nrow = 12))[1:5] / calib.comp$white_sim[calib.comp$par == "ODdeaths"]
  calib.comp$black_sim[calib.comp$par == "Fx_ODD"] <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "black"], nrow = 12))[1:5] / calib.comp$black_sim[calib.comp$par == "ODdeaths"]
  calib.comp$hisp_sim[calib.comp$par == "Fx_ODD"]  <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "hisp"], nrow = 12))[1:5] / calib.comp$hisp_sim[calib.comp$par == "ODdeaths"]
  calib.comp$white_sim[calib.comp$par == "EDvisits"] <- colSums(matrix(sim_sq$m.EDvisits.race[, "white"], nrow = 12))[1:5]
  calib.comp$black_sim[calib.comp$par == "EDvisits"] <- colSums(matrix(sim_sq$m.EDvisits.race[, "black"], nrow = 12))[1:5]
  calib.comp$hisp_sim[calib.comp$par == "EDvisits"]  <- colSums(matrix(sim_sq$m.EDvisits.race[, "hisp"], nrow = 12))[1:5]
  # 
  # write.xlsx(calib.comp, file = paste0("Outputs/MA/CalibrationPreliminary_", ss,".xlsx"), row.names = F) 
  # 
  # # mx.costs.totl[ss, "Status Quo"] <- sim_sq$total.cost
  # array.od.death.rgn[ , "white", "Status Quo"]    <- colSums(sim_sq$m.oddeath.white[(timesteps-23):timesteps, ])/2
  # array.od.death.rgn[ , "black", "Status Quo"]    <- colSums(sim_sq$m.oddeath.black[(timesteps-23):timesteps, ])/2
  # array.od.death.rgn[ , "hispanic", "Status Quo"] <- colSums(sim_sq$m.oddeath.hisp[(timesteps-23):timesteps, ])/2
  # array.oend.nlx.rgn[ , , "Status Quo"] <- sim_sq$n.nlx.OEND.str
  
  for (jj in 2:length(scenario.name)){
    params$Nx2ODratio <- ratio.level[jj-1]
    sim_pg <- MicroSim(init_ppl, params = params, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = scenario.name[jj], seed = seed) # run for program scenario
    result.cost[ss, jj] <- sim_pg$total.cost
    result.qaly[ss, jj] <- sim_pg$total.qaly
    result.qaly.w[ss, jj] <- sim_pg$total.qaly.w
    result.qaly.b[ss, jj] <- sim_pg$total.qaly.b
    result.qaly.h[ss, jj] <- sim_pg$total.qaly.h
    # mx.od.death.1st["total", jj]  <- sum(sim_pg$m.oddeath[(timesteps-35):(timesteps-24), ])
    # mx.od.death.1st["white", jj]  <- sum(sim_pg$m.oddeath.white[(timesteps-35):(timesteps-24), ])
    # mx.od.death.1st["black", jj]  <- sum(sim_pg$m.oddeath.black[(timesteps-35):(timesteps-24), ])
    # mx.od.death.1st["hisp", jj]   <- sum(sim_pg$m.oddeath.hisp[(timesteps-35):(timesteps-24), ])
    # mx.od.death.2nd["total", jj]  <- sum(sim_pg$m.oddeath[(timesteps-23):(timesteps-12), ])
    # mx.od.death.2nd["white", jj]  <- sum(sim_pg$m.oddeath.white[(timesteps-23):(timesteps-12), ])
    # mx.od.death.2nd["black", jj]  <- sum(sim_pg$m.oddeath.black[(timesteps-23):(timesteps-12), ])
    # mx.od.death.2nd["hisp", jj]   <- sum(sim_pg$m.oddeath.hisp[(timesteps-23):(timesteps-12), ])
    # mx.od.death.last["total", jj] <- sum(sim_pg$m.oddeath[(timesteps-11):timesteps, ])
    # mx.od.death.last["white", jj] <- sum(sim_pg$m.oddeath.white[(timesteps-11):timesteps, ])
    # mx.od.death.last["black", jj] <- sum(sim_pg$m.oddeath.black[(timesteps-11):timesteps, ])
    # mx.od.death.last["hisp", jj]  <- sum(sim_pg$m.oddeath.hisp[(timesteps-11):timesteps, ])
    # mx.od.death.totl["total", jj] <- sum(sim_pg$m.oddeath[(timesteps-35):timesteps, ])
    # mx.od.death.totl["white", jj] <- sum(sim_pg$m.oddeath.white[(timesteps-35):timesteps, ])
    # mx.od.death.totl["black", jj] <- sum(sim_pg$m.oddeath.black[(timesteps-35):timesteps, ])
    # mx.od.death.totl["hisp", jj]  <- sum(sim_pg$m.oddeath.hisp[(timesteps-35):timesteps, ])
    # mx.od.death.wtns.1st["total", jj]  <- sum(sim_pg$v.oddeath.w[(timesteps-35):(timesteps-24)])
    # mx.od.death.wtns.1st["white", jj]  <- sum(sim_pg$v.oddeath.w.race[(timesteps-35):(timesteps-24), "white"])
    # mx.od.death.wtns.1st["black", jj]  <- sum(sim_pg$v.oddeath.w.race[(timesteps-35):(timesteps-24), "black"])
    # mx.od.death.wtns.1st["hisp", jj]   <- sum(sim_pg$v.oddeath.w.race[(timesteps-35):(timesteps-24), "hisp"])
    # mx.od.death.wtns.2nd["total", jj]  <- sum(sim_pg$v.oddeath.w[(timesteps-23):(timesteps-12)])
    # mx.od.death.wtns.2nd["white", jj]  <- sum(sim_pg$v.oddeath.w.race[(timesteps-23):(timesteps-12), "white"])
    # mx.od.death.wtns.2nd["black", jj]  <- sum(sim_pg$v.oddeath.w.race[(timesteps-23):(timesteps-12), "black"])
    # mx.od.death.wtns.2nd["hisp", jj]   <- sum(sim_pg$v.oddeath.w.race[(timesteps-23):(timesteps-12), "hisp"])
    # mx.od.death.wtns.last["total", jj] <- sum(sim_pg$v.oddeath.w[(timesteps-11):timesteps])
    # mx.od.death.wtns.last["white", jj] <- sum(sim_pg$v.oddeath.w.race[(timesteps-11):timesteps, "white"])
    # mx.od.death.wtns.last["black", jj] <- sum(sim_pg$v.oddeath.w.race[(timesteps-11):timesteps, "black"])
    # mx.od.death.wtns.last["hisp", jj]  <- sum(sim_pg$v.oddeath.w.race[(timesteps-11):timesteps, "hisp"])
    # mx.od.death.wtns.totl["total", jj] <- sum(sim_pg$v.oddeath.w[(timesteps-35):timesteps])
    # mx.od.death.wtns.totl["white", jj] <- sum(sim_pg$v.oddeath.w.race[(timesteps-35):timesteps, "white"])
    # mx.od.death.wtns.totl["black", jj] <- sum(sim_pg$v.oddeath.w.race[(timesteps-35):timesteps, "black"])
    # mx.od.death.wtns.totl["hisp", jj]  <- sum(sim_pg$v.oddeath.w.race[(timesteps-35):timesteps, "hisp"])
    # # mx.costs.totl[ss, jj] <- sim_pg$total.cost
    # array.od.death.rgn[ , "white", jj]    <- colSums(sim_pg$m.oddeath.white[(timesteps-23):timesteps, ])/2
    # array.od.death.rgn[ , "black", jj]    <- colSums(sim_pg$m.oddeath.black[(timesteps-23):timesteps, ])/2
    # array.od.death.rgn[ , "hispanic", jj] <- colSums(sim_pg$m.oddeath.hisp[(timesteps-23):timesteps, ])/2
    # array.oend.nlx.rgn[ , , jj] <- sim_pg$n.nlx.OEND.str
  }
  
  # n.scenario <- length(scenario.name)
  # scenario <- c("Non-targeted", "Equity-focused")
  # sub.scenario <- c("40 kits per OOD", "60 kits per OOD", "80 kits per OOD")
  # rgn.results <- data.frame(location = rep(v.region, (n.scenario)), scenario = c(rep("Status Quo", num_regions), rep(scenario, each=num_regions*length(sub.scenario))), 
  #                           sub.scenario = c(rep("Status Quo", num_regions), rep(rep(sub.scenario, each = num_regions),length(scenario))),
  #                           N_Naloxone_total = numeric(num_regions * (n.scenario)), N_Naloxone_white = numeric(num_regions * (n.scenario)),
  #                           N_Naloxone_black = numeric(num_regions * (n.scenario)), N_Naloxone_hispanic = numeric(num_regions * (n.scenario)),
  #                           N_ODDeath_total = numeric(num_regions * (n.scenario)), N_ODDeath_white = numeric(num_regions * (n.scenario)),
  #                           N_ODDeath_black = numeric(num_regions * (n.scenario)), N_ODDeath_hispanic = numeric(num_regions * (n.scenario)),
  #                           population_total = rep(colSums(Demographic[(Demographic$race!="other"),-c(1:3)]), n.scenario),
  #                           population_white = rep(colSums(Demographic[(Demographic$race=="white"),-c(1:3)]), n.scenario),
  #                           population_black = rep(colSums(Demographic[(Demographic$race=="black"),-c(1:3)]), n.scenario),
  #                           population_hispanic = rep(colSums(Demographic[(Demographic$race=="hisp"),-c(1:3)]), n.scenario))
  # 
  # rgn.results$N_Naloxone_total[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- rowSums(array.oend.nlx.rgn[ , , "Status Quo"])
  # rgn.results$N_Naloxone_white[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- array.oend.nlx.rgn[ , "white", "Status Quo"]
  # rgn.results$N_Naloxone_black[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- array.oend.nlx.rgn[ , "black", "Status Quo"]
  # rgn.results$N_Naloxone_hispanic[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- array.oend.nlx.rgn[ , "hispanic", "Status Quo"]
  # rgn.results$N_ODDeath_total[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- rowSums(array.od.death.rgn[ , , "Status Quo"])
  # rgn.results$N_ODDeath_white[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- array.od.death.rgn[ , "white", "Status Quo"]
  # rgn.results$N_ODDeath_black[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- array.od.death.rgn[ , "black", "Status Quo"]
  # rgn.results$N_ODDeath_hispanic[rgn.results$scenario== "Status Quo" & rgn.results$sub.scenario == "Status Quo"] <- array.od.death.rgn[ , "hispanic", "Status Quo"]
  # 
  # for (ii in 1:length(scenario)){
  #   for (jj in 1:length(sub.scenario)){
  #     rgn.results$N_Naloxone_total[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- rowSums(array.oend.nlx.rgn[ , , (ii-1)*length(sub.scenario)+jj+1])
  #     rgn.results$N_Naloxone_white[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- array.oend.nlx.rgn[ , "white", (ii-1)*length(sub.scenario)+jj+1]
  #     rgn.results$N_Naloxone_black[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- array.oend.nlx.rgn[ , "black", (ii-1)*length(sub.scenario)+jj+1]
  #     rgn.results$N_Naloxone_hispanic[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- array.oend.nlx.rgn[ , "hispanic", (ii-1)*length(sub.scenario)+jj+1]
  #     rgn.results$N_ODDeath_total[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- rowSums(array.od.death.rgn[ , , (ii-1)*length(sub.scenario)+jj+1])
  #     rgn.results$N_ODDeath_white[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- array.od.death.rgn[ , "white", (ii-1)*length(sub.scenario)+jj+1]
  #     rgn.results$N_ODDeath_black[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- array.od.death.rgn[ , "black", (ii-1)*length(sub.scenario)+jj+1]
  #     rgn.results$N_ODDeath_hispanic[rgn.results$scenario== scenario[ii] & rgn.results$sub.scenario == sub.scenario[jj]] <- array.od.death.rgn[ , "hispanic", (ii-1)*length(sub.scenario)+jj+1]
  #   }
  # }
  # write.xlsx(mx.od.death.1st, file = paste0("Outputs/MA/ODdeaths_", ss,".xlsx"), sheetName = "1st year", row.names = T)
  # write.xlsx(mx.od.death.2nd, file = paste0("Outputs/MA/ODdeaths_", ss,".xlsx"), sheetName = "2nd year", append = T, row.names = T)
  # write.xlsx(mx.od.death.last, file = paste0("Outputs/MA/ODdeaths_", ss,".xlsx"), sheetName= "last year", append = T, row.names = T)
  # write.xlsx(mx.od.death.totl, file = paste0("Outputs/MA/ODdeaths_", ss,".xlsx"), sheetName= "Total", append = T, row.names = T)
  # write.xlsx(mx.od.death.wtns.1st, file = paste0("Outputs/MA/WitnessedDeaths_", ss,".xlsx"), sheetName = "1st year", row.names = T)
  # write.xlsx(mx.od.death.wtns.2nd, file = paste0("Outputs/MA/WitnessedDeaths_", ss,".xlsx"), sheetName = "2nd year", append = T, row.names = T)
  # write.xlsx(mx.od.death.wtns.last, file = paste0("Outputs/MA/WitnessedDeaths_", ss,".xlsx"), sheetName= "last year", append = T, row.names = T)
  # write.xlsx(mx.od.death.wtns.totl, file = paste0("Outputs/MA/WitnessedDeaths_", ss,".xlsx"), sheetName= "Total", append = T, row.names = T)
  # # write.xlsx(mx.costs.totl, file = ("Outputs/MA/Total_costs.xlsx"), row.names = F)
  # write.csv(rgn.results, file = paste0("Outputs/MA/Results_bytown_2y_", ss,".csv"), row.names = F)
}
write.csv(result.cost, file = ("Outputs/MA/Total_costs.csv"), row.names = F)
write.csv(result.qaly, file = ("Outputs/MA/Total_qalys.csv"), row.names = F)
write.csv(result.qaly.w, file = ("Outputs/MA/Total_qalys_white.csv"), row.names = F)
write.csv(result.qaly.b, file = ("Outputs/MA/Total_qalys_black.csv"), row.names = F)
write.csv(result.qaly.h, file = ("Outputs/MA/Total_qalys_hisp.csv"), row.names = F)

























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

n.scenario <- length(scenario.name)
rgn.results <- data.frame(location = rep(v.region, n.scenario), scenario = rep(scenario.name, each=num_regions), 
                          N_Nx_mean = numeric(num_regions * n.scenario), N_Nx_lower = numeric(num_regions * n.scenario), N_Nx_upper = numeric(num_regions * n.scenario),
                          N_ODDeath_mean = numeric(num_regions * n.scenario), N_ODDeath_lower = numeric(num_regions * n.scenario), N_ODDeath_upper = numeric(num_regions * n.scenario),
                          population = rep(colSums(Demographic[,-c(1:3)]), n.scenario))

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
