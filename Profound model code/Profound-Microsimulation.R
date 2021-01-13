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
rm(list=ls())
# install.packages("rstudioapi")
library(rstudioapi)
library(dplyr)
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
source("Profound-Function-PopulationCreation.R")
source("Profound-Function-TransitionProbability.R")
source("Profound-DecisionTree.R")

# INPUT PARAMETERS
set.seed(2020)                          # set the seed 

# Model structure
v.state     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
n.state     <- length(v.state)                                                     # number of states
n.t         <- 60                                                                  # number of time cycles (in month)
n.pop.t     <- 1000000                                                             # total population size for a given region
v.od        <- rep(0, times = n.t)                                                 # count of overdose events at each time step
v.oddeath   <- rep(0, times = n.t)                                                 # count of overdose deaths at each time step
v.str       <- c("SQ", "Expand50", "All+200", "R2+600")                            # store the strategy names
v.rgn       <- c("R1", "R2", "R3")                                                 # store the region names
n.rgn       <- length(v.rgn)                                                       # number of regions
d.c         <- 0.03                                                                # discounting of costs by 3%
cost.item   <- c("total", "naloxone")
cost.matrix <- matrix(0, nrow=n.t, ncol = length(cost.item))
colnames(cost.matrix) <- cost.item


# Parameters for initial cohort
prop.pwud      <- 0.008                  
mean.age       <- 40
prop.preb      <- 0.5
prop.illicit   <- 0.3
prop.il.lr     <- 0.7
prop.il.hr     <- 0.3
prop.inact     <- 0.1
prop.NOUD      <- 0.1
prop.fx        <- 0.2
prop.rsd       <- t(c(0.3, 0.6, 0.1));   colnames(prop.rsd) <- v.rgn
pop.prop.list  <- list(prop.illicit = prop.illicit, prop.il.lr = prop.il.lr, prop.il.hr = prop.il.hr,
                       prop.preb    = prop.preb   , prop.inact = prop.inact, prop.NOUD  = prop.NOUD,
                       prop.fx      = prop.fx     , prop.rsd   = prop.rsd)
n.pwud         <- n.pop.t * prop.pwud


# Initialize the study population - PWUD
pop.info  <- c("age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
init.pop  <- data.frame(matrix(0, n.pwud, length(pop.info)+1))                 # initial population matrix
dimnames(init.pop)[[2]] <- c("ind", pop.info)
init.pop$ind            <- 1:n.pwud
init.pop                <- pop.creat(init.pop, mean.age = mean.age, pop.prop.list, seed=2020)


# Overdose probability matrix (per month)
od.matrix             <- matrix(0, nrow = 4, ncol = 4)
rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
colnames(od.matrix)   <- c("first", "subs", "relap", "fx")
od.matrix["preb", "first"]   <- 0.001
od.matrix["il.lr", "first"]  <- 0.003
od.matrix["il.hr", "first"]  <- 0.012
od.matrix["NODU", "first"]   <- 0.001
multi.od2                    <- 2
multi.relap                  <- 3
multi.fx                     <- 5
od.matrix[ , "subs"]         <- od.matrix[ , "first"] * multi.od2
od.matrix[ , "relap"]        <- od.matrix[ , "first"] * multi.relap
od.matrix[ , "fx"]           <- od.matrix[ , "first"] * multi.fx


# Baseline mortality excluding overdose (per month)
mor.bl1.ut  <- 0.0187/12      # Age < 35 (ag1), not treated with MOUD
mor.bl1.t   <- 0.0073/12      # Age < 35 (ag1), treated with MOUD (inactive)
mor.bl2.ut  <- 0.0298/12      # Age > 35 (ag2), not treated with MOUD
mor.bl2.t   <- 0.0113/12      # Age > 35 (ag2), treated with MOUD (inactive)
mor.matrix                     <- matrix(0, nrow = 2, ncol = 2)
rownames(mor.matrix)           <- c("untreated", "treated")
colnames(mor.matrix)           <- c("ag1", "ag2")
mor.matrix["untreated", "ag1"] <- mor.bl1.ut
mor.matrix["treated", "ag1"]   <- mor.bl1.t
mor.matrix["untreated", "ag2"] <- mor.bl2.ut
mor.matrix["treated", "ag2"]   <- mor.bl2.t


# Transistion probabilities between states (per month)
p.preb2il.lr  <- 0.0034             # probability from prescribed to illicit, low-risk
p.preb2inact  <- 0.00903            # probability from prescribed to inactive
p.il.lr2il.hr <- 0.0049          	  # probability from illicit, low-risk to high-risk
p.il.lr2inact <- 0.0051          	  # probability from illicit, low-risk to inactive
p.il.hr2il.lr <- 0.0025          	  # probability from illicit, high-risk to low-risk
p.il.hr2inact <- 0.0051          	  # probability from illicit, high-risk to inactive
p.inact2relap <- 0.07               # probability to replased opioid use


# Naloxone distribution
n.nlx.v       <- c(0, 1000, 200)
n.od_death.v  <- c(150, 800, 250)
nlx.adj       <- 1                  # adjuster for naloxone availability given the ratio between naloxone kits to od_deaths in a region

# Spatial dynamics to access naloxone kits
R1.nlx  <-  c(0.7, 0.2, 0.1)
R2.nlx  <-  c(0.15, 0.8, 0.05)
R3.nlx  <-  c(0.1, 0.25, 0.65)
acc.nlx.matrix <- rbind(R1.nlx, R2.nlx, R3.nlx)
rownames(acc.nlx.matrix)  <- c("R1.from", "R2.from", "R3.from")
colnames(acc.nlx.matrix)  <- c("R1.to", "R2.to", "R3.to")


# Cost inputs for each health state
c.preb             <- 800                    # cost of remaining one cycle: prescribed OU
c.il.lr            <- 1019                   # cost of remaining one cycle: illicit OU, low-risk
c.il.hr            <- 1511                   # cost of remaining one cycle: illicit OU, high-risk
c.inact            <- 631                    # cost of remaining one cycle: inactive OU
c.NODU             <- 2000                   # cost of remaining one cycle: non-opioid
c.relap.v          <- numeric(0)             # cost of remaining one cycle: relapsed, as the average of inactive and prior state
c.relap.v["preb"]  <- (c.preb + c.inact)/2
c.relap.v["il.lr"] <- (c.il.lr + c.inact)/2
c.relap.v["il.hr"] <- (c.il.hr + c.inact)/2
c.nlx.kit          <- 70                     # unit cost for naloxone per kit
c.nlx.dtb          <- 10                     # unit service cost per naloxone kit distributed

# Cost inputs for decision tree nodes
c.EMS              <- 2092                   # cost for EMS visit
c.hospcare         <- 1034                   # cost of hospital care for overdosed PWUD admitted to hospital


MicroSim <- function(init.pop, n.pwud, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = 1) {
  # Arguments:  
  # init.pop:      matrix of initial states for individuals
  # n.pwud:        number of PWUD
  # n.t:           total number of cycles to run the model
  # v.state:       vector of health state names
  # d.c:           discount rate for costs
  # PT.out:        should the output include a Microsimulation trace? (default is TRUE)
  # Str:           simulating strategy
  # seed:          starting seed number for random number generator (default is 1)
  # Makes use of:
  # trans.prob:    function for the estimation of transition probabilities
  # trans.prob:    function for the estimation of transition probabilities
  # Costs:         function for the estimation of cost state values
  # decision.tree: function for the decision tree module  
  
  if (Str == "SQ"){
    n.nlx.v <- n.nlx.v
  } else if (Str == "Expand50"){
    n.nlx.v <- n.nlx.v + n.nlx.v * 0.5
  } else if (Str == "All+200"){
    n.nlx.v <- n.nlx.v + 200
  } else if (Str == "R2+600"){
    n.nlx.v[2] <- n.nlx.v[2] + 600
  }
  
  v.dwc <- 1 / (1 + d.c) ^ (0:(n.t-1))   # calculate the cost discount weight based on the discount rate d.c
  
  # Create the population list to capture the state/attributes/costs for all individuals at each time point 
  pop.list <- list()
  pop.list[[1]] <- init.pop       # indicate the initial health state and attributes
  
  set.seed(seed)                  # set the seed for every individual for the random number generator

  cost.matrix[1, ] <- Costs(state = pop.list[[1]]$curr.state, OU.state = pop.list[[1]]$OU.state, nlx = sum(n.nlx.v)/12 , count = NULL)  # estimate costs per individual for the initial health state]
  
  for (t in 2:n.t) {
    m.tp                       <- trans.prob(pop.list[[t-1]])                # calculate the transition probabilities at cycle t 
    pop.list[[t]]              <- pop.list[[t-1]]
    pop.list[[t]]$curr.state   <- as.vector(samplev(probs = m.tp, m = 1))  # sample the next health state and store that state in matrix m.M
    ind.oustate.chg            <- filter(pop.list[[t]], curr.state %in% v.oustate &  OU.state != curr.state)$ind
    pop.list[[t]]$OU.state[ind.oustate.chg] <- pop.list[[t]]$curr.state[ind.oustate.chg]
    
    od.pop                     <- pop.list[[t]][pop.list[[t]]$curr.state == "od", ]
    v.od[t]                    <- nrow(od.pop)
    
    decntree.out               <- decision.tree(od.pop, n.nlx.v, n.od_death.v, acc.nlx.matrix, seed = seed)
    
    v.oddeath[t]               <- sum(decntree.out[ , "od.death"])
    n.EMS                      <- sum(decntree.out[ , "EMS"])
    n.hospcare                 <- sum(decntree.out[ , "hospcare"])
    od.pop$curr.state[decntree.out[ , "od.death"] == 1]                           <- "dead"
    od.pop$ever.od[decntree.out[ , "od.death"] != 1]                              <- 1
    od.pop$curr.state[decntree.out[ , "inact"] == 1]                              <- "inact"
    od.pop$curr.state[od.pop$curr.state != "dead" | od.pop$curr.state != "dead"]  <- od.pop$OU.state[od.pop$curr.state != "dead" | od.pop$curr.state != "dead"]
    
    pop.list[[t]][od.pop$ind, ] <- od.pop
    cost.matrix[t, ]  <- Costs(state = pop.list[[t]]$curr.state, OU.state = pop.list[[t]]$OU.state, nlx = sum(n.nlx.v)/12 , count = list(n.EMS = n.EMS, n.hospcare = n.hospcare))

    pop.list[[t]]$age <- pop.list[[t]]$age + 1 
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
  } # close the loop for the time points 
  
  
  total.cost <- sum(cost.matrix[ , "total"] * v.dwc)       # total (discounted) cost
  
  if (PT.out == TRUE){
    pop.trace = pop.list
  } else{
    pop.trace = NULL
  }
  
  results <- list(v.oddeath = v.oddeath, v.od = v.od, cost.matrix = cost.matrix, total.cost = total.cost, pop.trace = pop.trace) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  


### Costs function
# The Costs function estimates the costs at every cycle.

Costs <- function (state, OU.state, nlx, count) {
  if (is.null(count)){
    count.EMS       <- 0
    count.hospcare  <- 0
  } else{
    count.EMS       <- count$n.EMS
    count.hospcare  <- count$n.hospcare
  }
  c.TC  <- sum(state == "preb")  * c.preb     +
           sum(state == "il.lr") * c.il.lr    +
           sum(state == "il.hr") * c.il.hr    +
           sum(state == "inact") * c.inact    +
           sum(state == "NODU")  * c.NODU     +
           sum(state == "relap" & OU.state == "preb")  * c.relap.v["preb"]  +
           sum(state == "relap" & OU.state == "il.lr") * c.relap.v["il.lr"]  +
           sum(state == "relap" & OU.state == "il.hr") * c.relap.v["il.hr"]  +
           count.EMS       * c.EMS            +
           count.hospcare  * c.hospcare       +
           nlx * (c.nlx.dtb + c.nlx.kit)
           
  c.nlx  <- nlx * (c.nlx.dtb + c.nlx.kit)
  return(c(c.TC, c.nlx))        		                   # return the costs
}



##################################### Run the simulation ##################################
# START SIMULATION
p = Sys.time()
sim_sq   <- MicroSim(init.pop, n.pwud, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = 100) # run for no treatment
sim_ep50 <- MicroSim(init.pop, n.pwud, n.t, v.state, d.c, PT.out = TRUE, Str = "Expand50", seed = 100) # run for treatment
sim_a200 <- MicroSim(init.pop, n.pwud, n.t, v.state, d.c, PT.out = TRUE, Str = "All+200", seed = 100) # run for treatment
sim_r2   <- MicroSim(init.pop, n.pwud, n.t, v.state, d.c, PT.out = TRUE, Str = "R2+600", seed = 100) # run for treatment

comp.time = Sys.time() - p

comp.time

sum(sim_sq$v.oddeath)
sum(sim_ep50$v.oddeath)
sum(sim_a200$v.oddeath)
sum(sim_r2$v.oddeath)