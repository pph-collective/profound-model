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


MicroSim <- function(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = 1) {
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
  # Costs:         function for the estimation of cost state values
  # decision.tree: function for the decision tree module  
  
  n.nlx.mx.lst  <- array.Nx[dim(array.Nx)[1], , ]
  if (Str == "SQ"){
    n.nlx.mx.str <- n.nlx.mx.lst
  } else if (Str == "Expand50"){
    n.nlx.mx.str <- n.nlx.mx.lst + n.nlx.mx.lst * 0.5
  } else if (Str == "All+200"){
    n.nlx.mx.str <- n.nlx.mx.lst + 200
  } else if (Str == "R2+600"){
    n.nlx.mx.str[2] <- n.nlx.mx.lst[2] + 600
  }
  array.Nx <- abind(array.Nx, n.nlx.mx.str, along = 1)
  
  v.dwc <- rep(1 / (1 + d.c) ^ (0:(n.yr-1)), each =12)   # calculate the cost discount weight based on the discount rate d.c
  
  # Create the population list to capture the state/attributes/costs for all individuals at each time point 
  pop.list <- list()
  # pop.list[[1]] <- init.pop       # indicate the initial health state and attributes
  
  set.seed(seed)                  # set the seed for every individual for the random number generator
  
  # cost.matrix[1, ] <- Costs(state = pop.list[[1]]$curr.state, OU.state = pop.list[[1]]$OU.state, nlx = sum(array.Nx[1,,])/12 , count = NULL)  # estimate costs per individual for the initial health state]
  
  for (t in 1:n.t) {
    n.nlx.yr                     <- array.Nx[floor((t-1)/12)+1, , ]
    if (t == 1){
      m.tp                       <- trans.prob(init.pop)                # calculate the transition probabilities at cycle t 
      pop.list[[t]]              <- init.pop
      n.nlx.mn                   <- init.Nx + n.nlx.yr/12
    } else {
      m.tp                       <- trans.prob(pop.list[[t-1]])                # calculate the transition probabilities at cycle t 
      pop.list[[t]]              <- pop.list[[t-1]]
      n.nlx.mn                   <- n.nlx.mn*(1-r.LossExp) + n.nlx.yr/12
    }
    pop.list[[t]]$curr.state   <- as.vector(samplev(probs = m.tp, m = 1))    # sample the next health state and store that state in matrix m.M
    ind.oustate.chg            <- filter(pop.list[[t]], curr.state %in% v.oustate &  OU.state != curr.state)$ind
    pop.list[[t]]$OU.state[ind.oustate.chg] <- pop.list[[t]]$curr.state[ind.oustate.chg]
    
    od.pop                     <- pop.list[[t]][pop.list[[t]]$curr.state == "od", ]
    v.od[t]                    <- nrow(od.pop)
    ou.pop.resid               <- pop.list[[t]] %>% count(residence)
    
    decntree.out               <- decision.tree(od.pop, n.nlx = n.nlx.mn, ou.pop.resid, seed = seed+t)
    
    v.oddeath[t]               <- sum(decntree.out[ , "od.death"])
    n.EMS                      <- sum(decntree.out[ , "EMS"])
    n.hospcare                 <- sum(decntree.out[ , "hospcare"])
    od.pop$curr.state[decntree.out[ , "od.death"] == 1]   <- "dead"
    od.pop$ever.od[decntree.out[ , "od.death"] != 1]      <- 1
    od.pop$curr.state[decntree.out[ , "inact"] == 1]      <- "inact"
    od.pop$curr.state[od.pop$curr.state != "dead"]        <- od.pop$OU.state[od.pop$curr.state != "dead"]
    od.death.sum <- od.pop[od.pop$curr.state == "dead", ] %>% count(residence)
    for (dd in 1:nrow(od.death.sum)){
      m.oddeath[t , od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
    pop.list[[t]][od.pop$ind, ] <- od.pop
    cost.matrix[t, ]  <- Costs(state = pop.list[[t]]$curr.state, OU.state = pop.list[[t]]$OU.state, nlx = sum(n.nlx.yr)/12 , count = list(n.EMS = n.EMS, n.hospcare = n.hospcare))
    
    pop.list[[t]]$age[pop.list[[t]]$curr.state != "dead"] <- pop.list[[t]]$init.age[pop.list[[t]]$curr.state != "dead"] + floor(t/12) 
    cat('\r', paste(round(t/n.t * 100), "% done", sep = " "))       # display the progress of the simulation
  } # close the loop for the time points 
  
  
  total.cost <- sum(cost.matrix[ , "TotalCost"] * v.dwc)       # total (discounted) cost
  
  if (PT.out == TRUE){
    pop.trace = pop.list
  } else{
    pop.trace = NULL
  }
  
  results <- list(v.oddeath = v.oddeath, m.oddeath = m.oddeath, v.od = v.od, cost.matrix = cost.matrix, total.cost = total.cost, pop.trace = pop.trace) # store the results from the simulation in a list  
  return(results)  # return the results
}  # end of the MicroSim function  
