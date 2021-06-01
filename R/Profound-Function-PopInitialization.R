########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Function module for the microsimulation of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: May 16, 2020
# Last update: May 30, 2020
#
########################################################################################
#################        Population creation function      #############################
########################################################################################

pop.initiation <- function(initials, seed = 2021){
  list2env(initials, environment())
  ## opioid population
  oud_pop_nsduh <- demo.mx * OUDDemo                              #number of OUD population estimates for each subgroup in each region according to NSDUH prevalence estimates
  oud.adj       <- prev.oud / (sum(oud_pop_nsduh)/sum(demo.mx))   #adjuster for OUD population according to overall oud prevalence and estimated one from the previous step
  oud_pop       <- oud_pop_nsduh * oud.adj
  prop.oud.rgn  <- colSums(oud_pop) / sum(oud_pop)
  oud_demo      <- oud_pop / matrix(rep(colSums(oud_pop), each = nrow(oud_pop)), nrow = nrow(oud_pop), ncol = ncol(oud_pop))
  
  ## non-opioid population
  stim_pop      <- demo.mx * StimDemo                              #number of NODU population estimates for each subgroup in each region according to NSDUH prevalence estimates
  prop.nodu.rgn <- colSums(stim_pop) / sum(stim_pop)
  nodu_demo     <- stim_pop / matrix(rep(colSums(stim_pop), each = nrow(stim_pop)), nrow = nrow(stim_pop), ncol = ncol(stim_pop))
  
  ## initialize the population matrix
  total.oud     <- round(sum(oud_pop), 0)
  total.nodu    <- round(sum(stim_pop), 0)
  init.pop      <- data.frame(matrix(0, (total.oud + total.nodu), length(pop.info)+1))                 # initial population matrix
  dimnames(init.pop)[[2]] <- c("ind", pop.info)
  init.pop$ind            <- 1:(total.oud + total.nodu)
  
  set.seed(seed)
  oud.rgn.ind   <- sample(1:length(v.rgn), size = total.oud,  prob = prop.oud.rgn,  replace = TRUE)
  nodu.rgn.ind  <- sample(1:length(v.rgn), size = total.nodu, prob = prop.nodu.rgn, replace = TRUE)
  
  ## OPIOID USE POPULATION
  for (i in 1:total.oud){
    set.seed(seed+i)
    # determine demographic
    demo.ind   <- sample(1:nrow(oud_demo), size = 1,  prob = oud_demo[,oud.rgn.ind[i]])
    sex        <- v.demo.sex[demo.ind]
    race       <- v.demo.race[demo.ind]
    age.ind    <- v.demo.age[demo.ind]
    if(age.ind == "12to17"){
      age <- init.age <- sample(12:17, size = 1)
    } else if(age.ind == "18to25"){
      age <- init.age <- sample(18:25, size = 1)
    } else if(age.ind == "26to34"){
      age <- init.age <- sample(26:34, size = 1)
    } else if(age.ind == "35to49"){
      age <- init.age <- sample(35:49, size = 1)
    } else if(age.ind == "50to64"){
      age <- init.age <- sample(50:64, size = 1)
    } else if(age.ind == "65older"){
      age <- init.age <- 65
    }
    
    # determine opioid use state
    oud.state  <- v.state[1:4]
    oud.prob.m <- c((1-ini.inact)* (1-ini.il.m), (1-ini.inact)* ini.il.m* (1-ini.il.hr.m), (1-ini.inact)* ini.il.m* ini.il.hr.m, ini.inact)
    oud.prob.f <- c((1-ini.inact)* (1-ini.il.f), (1-ini.inact)* ini.il.f* (1-ini.il.hr.f), (1-ini.inact)* ini.il.f* ini.il.hr.f, ini.inact)
    if (sex == "m"){
      curr.state <- init.state <- oud.state[sample(1:length(oud.state), size = 1, prob = oud.prob.m)]
    } else {
      curr.state <- init.state <- oud.state[sample(1:length(oud.state), size = 1, prob = oud.prob.f)]
    }

    if (curr.state == "inact"){
      if (sex == "m"){
        OU.state <- oud.state[sample(1:length(oud.state[1:3]), size = 1, prob = oud.prob.m[1:3]/sum(oud.prob.m[1:3]))]
      } else {
        OU.state <- oud.state[sample(1:length(oud.state[1:3]), size = 1, prob = oud.prob.f[1:3]/sum(oud.prob.f[1:3]))]
      }
    } else {
      OU.state <- curr.state
    }
    
    # # determine fentanyl use
    # fx         <- sample(0:1, size = 1, prob = c(1-ini.OUD.fx, ini.OUD.fx))
    
    # determine ever overdosed
    if (OU.state == "preb"){
      ever.od    <- sample(0:1, size = 1, prob = c(1-ini.everod.preb, ini.everod.preb))
    } else if (OU.state == "il.lr"){
      ever.od    <- sample(0:1, size = 1, prob = c(1-ini.everod.il.lr, ini.everod.il.lr))
    } else if (OU.state == "il.hr"){
      ever.od    <- sample(0:1, size = 1, prob = c(1-ini.everod.il.hr, ini.everod.il.hr))
    }

    
    init.pop$sex[i]        <- sex
    init.pop$race[i]       <- race
    init.pop$age[i]        <- age
    init.pop$residence[i]  <- v.rgn[oud.rgn.ind[i]]
    init.pop$curr.state[i] <- curr.state
    init.pop$OU.state[i]   <- OU.state
    init.pop$init.age[i]   <- init.age
    init.pop$init.state[i] <- init.state
    init.pop$ever.od[i]    <- ever.od
    # init.pop$fx[i]         <- fx
  }
  
  ## NON-OPIOID USE POPULATION
  for (i in 1:total.nodu){
    set.seed(seed+i)
    
    curr.state <- init.state <- OU.state <- "NODU"
    
    # determine demographic
    demo.ind   <- sample(1:nrow(nodu_demo), size = 1,  prob = nodu_demo[,nodu.rgn.ind[i]])
    sex        <- v.demo.sex[demo.ind]
    race       <- v.demo.race[demo.ind]
    age.ind    <- v.demo.age[demo.ind]
    if(age.ind == "12to17"){
      age <- init.age <- sample(12:17, size = 1)
    } else if(age.ind == "18to25"){
      age <- init.age <- sample(18:25, size = 1)
    } else if(age.ind == "26to34"){
      age <- init.age <- sample(26:34, size = 1)
    } else if(age.ind == "35to49"){
      age <- init.age <- sample(35:49, size = 1)
    } else if(age.ind == "50to64"){
      age <- init.age <- sample(50:64, size = 1)
    } else if(age.ind == "65older"){
      age <- init.age <- 65
    }
    
    # # determine fentanyl use
    # fx         <- sample(0:1, size = 1, prob = c(1-ini.NOUD.fx, ini.NOUD.fx))
    
    # determine ever overdosed
    ever.od    <- sample(0:1, size = 1, prob = c(1-ini.everod.sti, ini.everod.sti))
    
    
    init.pop$sex[i+total.oud]        <- sex
    init.pop$race[i+total.oud]       <- race
    init.pop$age[i+total.oud]        <- age
    init.pop$residence[i+total.oud]  <- v.rgn[nodu.rgn.ind[i]]
    init.pop$curr.state[i+total.oud] <- curr.state
    init.pop$OU.state[i+total.oud]   <- OU.state
    init.pop$init.age[i+total.oud]   <- init.age
    init.pop$init.state[i+total.oud] <- init.state
    init.pop$ever.od[i+total.oud]    <- ever.od
    # init.pop$fx[i+total.oud]         <- fx
  }
  return(init.pop)
}