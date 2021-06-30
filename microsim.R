###############################################################################################
#########################           Microsimulation        ####################################
###############################################################################################

###############################################################################################
####    Microsimulation to determine health states and number of overdoses                 ####
####    6 health states: prescribed, illicit (L/H), inactive, non-opioid, relapsed, death  ####
####    1 health event:  Overdose                                                          ####
####    Attributes:      state, age, sex, fentanyl, overdosed, pre.state,                  ####
####    Built to inform Naloxone distribution strategies to prevent overdsoe death         ####
###############################################################################################


MicroSim <- function(init.pop, vparameters, n.t, v.state, discount.rate, PT.out = TRUE, Str = "SQ", seed = 1) {
  # Arguments:  
  # init.pop:       matrix of initial states for individuals
  # vparameters:    ???
  # n.t:            total number of cycles to run the model
  # v.state:        vector of health state names
  # discount.rate:  discount rate for costs
  # PT.out:         should the output include a Microsimulation trace? (default is TRUE)
  # Str:            simulating strategy
  # seed:           starting seed number for random number generator (default is 1)
  # Makes use of:  TO_REVIEW: is this a normal thing to put in a docstring
  # trans.prob:     function for the estimation of transition probabilities
  # Costs:          function for the estimation of cost state values
  # decision.tree:  function for the decision tree module
  # TODO: actual docstring description
  list2env(vparameters, environment())
  # Find number of opioid and non-opioid users
  n.opioid <- sum(init.pop$curr.state != "NODU")
  n.noud   <- sum(init.pop$curr.state == "NODU")
  init.pop.residence   <- (init.pop %>% count(residence))$n
  print(NxDataPharm)
  # Create matrixes for ??? TO_REVIEW
  # TO_REVIEW: what is NxDataPharm? It's not clear from the name, and the vparameters make it difficult to track down
  NxPharm.mx      <- NxDataPharm$pe[NxDataPharm$year>=(yr.first-1)] %*% t(init.pop.residence / sum(init.pop.residence))
  NxPharm.array   <- array(0, dim = c(dim(NxPharm.mx)[1], 2, dim(NxPharm.mx)[2]))
  for (cc in 1:dim(NxPharm.mx)[1]){
    NxPharm.array[ cc, , ] <- round(rep(NxPharm.mx[cc, ], each = 2) * OD_loc, 0)
  }
  
  array.Nx    <- NxOEND.array[dimnames(NxOEND.array)[[1]] >= yr.first,   , ] + NxPharm.array[ -1, , ]
  init.Nx     <- NxOEND.array[dimnames(NxOEND.array)[[1]] == yr.first-1, , ] + NxPharm.array[  1, , ]
  
  n.nlx.mx.lst  <- array.Nx[dim(array.Nx)[1], , ]
  if (Str == "SQ"){
    n.nlx.mx.str <- n.nlx.mx.lst
  } else if (Str == "expand"){ 
    n.nlx.mx.str <- NxOEND.array[dim(NxOEND.array)[1],   , ] * 2 + NxPharm.array[dim(NxPharm.array)[1],   , ]
  } else if (Str == "program"){
    n.nlx.mx.str <- n.nlx.mx.lst + pg.add
  }
  
  array.Nx <- abind(array.Nx, n.nlx.mx.str, along = 1)
  
  v.dwc <- rep(1 / (1 + discount.rate) ^ (0:(n.yr-1)), each =12)   # calculate the cost discount weight based on the discount rate
  
  # Create the population list to capture the state/attributes/costs for all individuals at each time point 
  pop.list <- list()
  set.seed(seed)                  # set the seed for every individual for the random number generator
  
  for (t in 1:n.t) {
    n.nlx.yr                     <- array.Nx[floor((t-1)/12)+1, , ]
    if (t == 1){
      pop.list[[t]]              <- init.pop
      OUD.fx                     <- ini.OUD.fx
      # determine fentanyl use among initial population who use opioids 
      set.seed(seed)
      fx         <- sample(0:1, size = n.opioid, prob = c(1-OUD.fx, OUD.fx), replace = T)
      pop.list[[t]]$fx[init.pop$curr.state != "NODU"] <- fx
      # determine fentanyl use among initial population who use stimulants (non-opioid)
      set.seed(seed*2)
      fx         <- sample(0:1, size = n.noud, prob = c(1-ini.NOUD.fx, ini.NOUD.fx), replace = T)
      pop.list[[t]]$fx[init.pop$curr.state == "NODU"] <- fx
      m.tp                       <- trans.prob(pop.list[[t]], vparameters)                # calculate the transition probabilities at cycle t 
      n.nlx.mn                   <- init.Nx + n.nlx.yr/12
    } else {
      pop.list[[t]]              <- pop.list[[t-1]]
      if (t %% 12 ==0){
        OUD.fx     <- min(ini.OUD.fx * (1 + gw.fx * min(floor((t-1)/12)+1, 3)), 0.9)
        # determine fentanyl use among initial population who use opioids 
        set.seed(seed)
        fx         <- sample(0:1, size = n.opioid, prob = c(1-OUD.fx, OUD.fx), replace = T)
        pop.list[[t]]$fx[init.pop$curr.state != "NODU"] <- fx
        # # determine fentanyl exposure among population who use stimulants (non-opioid)
        # set.seed(seed*2)
        # fx         <- sample(0:1, size = n.noud, prob = c(1-ini.NOUD.fx, ini.NOUD.fx), replace = T)
        # pop.list[[t]]$fx[init.pop$curr.state == "NODU"] <- fx
      }
      m.tp                       <- trans.prob(pop.list[[t-1]], vparameters)                # calculate the transition probabilities at cycle t
      n.nlx.mn                   <- n.nlx.mn*(1-r.LossExp) + n.nlx.yr/12
    }
    pop.list[[t]]$curr.state   <- as.vector(samplev(probs = m.tp, m = 1))    # sample the next health state and store that state in matrix m.M
    ind.oustate.chg            <- filter(pop.list[[t]], curr.state %in% v.oustate &  OU.state != curr.state)$ind
    pop.list[[t]]$OU.state[ind.oustate.chg] <- pop.list[[t]]$curr.state[ind.oustate.chg]
    
    od.pop                     <- pop.list[[t]][pop.list[[t]]$curr.state == "od", ]
    v.od[t]                    <- nrow(od.pop)
    ou.pop.resid               <- pop.list[[t]] %>% count(residence)
    
    decntree.out               <- decision.tree(od.pop, n.nlx = n.nlx.mn, ou.pop.resid, vparameters, seed = seed+t)
    
    v.oddeath[t]               <- sum(decntree.out[ , "od.death"])
    v.odpriv[t]                <- sum(decntree.out[ , "locpriv"])
    v.odpubl[t]                <- v.od[t] - v.odpriv[t]
    v.deathpriv[t]             <- sum(decntree.out[decntree.out[ , "od.death"] == 1, "locpriv"])
    v.deathpubl[t]             <- sum(decntree.out[ , "od.death"] == 1) - v.deathpriv[t]
    v.nlxused[t]               <- sum(decntree.out[ , "nlx.used"])
    print("Here")
    n.EMS                      <- sum(decntree.out[ , "EMS"])
    n.hospcare                 <- sum(decntree.out[ , "hospcare"])
    od.pop$curr.state[decntree.out[ , "od.death"] == 1]   <- "dead"
    od.pop$ever.od[decntree.out[ , "od.death"] != 1]      <- 1
    od.pop$curr.state[decntree.out[ , "inact"] == 1]      <- "inact"
    od.pop$curr.state[od.pop$curr.state == "od"]          <- od.pop$OU.state[od.pop$curr.state == "od"]
    
    m.oddeath.fx[t] <- nrow(od.pop[od.pop$curr.state == "dead" & od.pop$fx ==1, ])
    m.oddeath.op[t] <- nrow(od.pop[od.pop$curr.state == "dead" & od.pop$OU.state != "NODU", ])
    m.oddeath.hr[t] <- nrow(od.pop[od.pop$curr.state == "dead" & od.pop$OU.state != "NODU" & od.pop$OU.state != "preb", ])
    m.oddeath.st[t] <- nrow(od.pop[od.pop$curr.state == "dead" & od.pop$OU.state == "NODU", ]) 
    m.EDvisits[t]   <- n.hospcare
      
    od.death.sum <- od.pop[od.pop$curr.state == "dead", ] %>% count(residence)
    for (dd in 1:nrow(od.death.sum)){
      m.oddeath[t , od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
    pop.list[[t]][od.pop$ind, ] <- od.pop
    cost.matrix[t, ]  <- Costs(state = pop.list[[t]]$curr.state, OU.state = pop.list[[t]]$OU.state, nlx = sum(n.nlx.yr)/12 , count = list(n.EMS = n.EMS, n.hospcare = n.hospcare), vparameters)
    
    pop.list[[t]]$age[pop.list[[t]]$curr.state != "dead"] <- pop.list[[t]]$init.age[pop.list[[t]]$curr.state != "dead"] + floor(t/12)     #update age for individuals that are still alive
    
    ##replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0)
    pop.list[[t]]$age[pop.list[[t]]$curr.state == "dead"]        <- pop.list[[t]]$init.age[pop.list[[t]]$curr.state == "dead"]
    pop.list[[t]]$ever.od[pop.list[[t]]$curr.state == "dead"]    <- 0
    pop.list[[t]]$OU.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state != "inact"]   <- pop.list[[t]]$init.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state != "inact"]
    pop.list[[t]]$OU.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state == "inact"]   <- pop.list[[t]]$OU.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state == "inact"]
    pop.list[[t]]$curr.state[pop.list[[t]]$curr.state == "dead"] <- pop.list[[t]]$init.state[pop.list[[t]]$curr.state == "dead"]
    
    # cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
  } # end the loop for the time steps 
  
  
  total.cost <- sum(cost.matrix[ , "TotalCost"] * v.dwc)       # total (discounted) cost
  
  if (PT.out == TRUE){
    pop.trace = pop.list
  } else{
    pop.trace = NULL
  }
  print("Saving results")
  results <- list(v.oddeath = v.oddeath, m.oddeath = m.oddeath, v.od = v.od, 
                  cost.matrix = cost.matrix, total.cost = total.cost, pop.trace = pop.trace, n.nlx.OEND.str = (n.nlx.mx.str - NxPharm.array[dim(NxPharm.array)[1],   , ]), n.nlx.all.str = n.nlx.mx.str,
                  m.oddeath.fx = m.oddeath.fx, m.oddeath.op = m.oddeath.op, m.oddeath.st = m.oddeath.st, m.oddeath.hr= m.oddeath.hr, m.EDvisits= m.EDvisits,
                  v.odpriv = v.odpriv, v.odpubl = v.odpubl, v.deathpriv = v.deathpriv, v.deathpubl = v.deathpubl, v.nlxused = v.nlxused) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  
