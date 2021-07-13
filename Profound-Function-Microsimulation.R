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


MicroSim <- function(init_ppl, params, timesteps, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = 1) {
  # Arguments:
  # init_ppl:      matrix of initial states for individuals
  # n.pwud:        number of PWUD
  # timesteps:           total number of cycles to run the model
  # v.state:       vector of health state names
  # d.c:           discount rate for costs
  # PT.out:        should the output include a Microsimulation trace? (default is TRUE)
  # Str:           simulating strategy
  # seed:          starting seed number for random number generator (default is 1)
  # Makes use of:
  # trans.prob:    function for the estimation of transition probabilities
  # Costs:         function for the estimation of cost state values
  # decisiotimesteptree: function for the decision tree module
  list2env(params, environment())
  n.opioid <- sum(init_ppl$curr.state != "NODU")
  n.noud <- sum(init_ppl$curr.state == "NODU")
  ini.pop.resid <- (init_ppl %>% count(residence))$n
  NxPharm.mx <- NxDataPharm$pe[NxDataPharm$year >= (yr.first - 1)] %*% t(ini.pop.resid / sum(ini.pop.resid))
  NxPharm.array <- array(0, dim = c(dim(NxPharm.mx)[1], 2, dim(NxPharm.mx)[2]))
  for (cc in 1:dim(NxPharm.mx)[1]) {
    NxPharm.array[cc, , ] <- round(rep(NxPharm.mx[cc, ], each = 2) * OD_loc, 0)
  }

  array.Nx <- NxOEND.array[dimnames(NxOEND.array)[[1]] >= yr.first, , ] + NxPharm.array[-1, , ]
  init.Nx <- NxOEND.array[dimnames(NxOEND.array)[[1]] == yr.first - 1, , ] + NxPharm.array[1, , ]

  n.nlx.mx.lst <- array.Nx[dim(array.Nx)[1], , ]
  if (Str == "SQ") {
    n.nlx.mx.str <- n.nlx.mx.lst
  } else if (Str == "expand") {
    n.nlx.mx.str <- NxOEND.array[dim(NxOEND.array)[1], , ] * exp.lv + NxPharm.array[dim(NxPharm.array)[1], , ]
  } else if (Str == "program") {
    n.nlx.mx.str <- n.nlx.mx.lst + pg.add
  }

  for (aa in 1:(n.yr - dim(array.Nx)[1])) {
    array.Nx <- abind(array.Nx, n.nlx.mx.str, along = 1)
  }

  v.dwc <- rep(1 / (1 + d.c)^(0:(n.yr - 1)), each = 12) # calculate the cost discount weight based on the discount rate d.c

  # Create the population list to capture the state/attributes/costs for all individuals at each time point
  pop.list <- list()
  set.seed(seed) # set the seed for every individual for the random number generator

  for (t in 1:timesteps) {
    n.nlx.yr <- array.Nx[floor((t - 1) / 12) + 1, , ]
    if (t == 1) {
      pop.list[[t]] <- init_ppl
      OUD.fx <- ini.OUD.fx
      # determine fentanyl use among population who use opioids
      set.seed(seed)
      fx <- sample(0:1, size = n.opioid, prob = c(1 - OUD.fx, OUD.fx), replace = T)
      pop.list[[t]]$fx[init_ppl$curr.state != "NODU"] <- fx
      # determine fentanyl use among population who use stimulants (non-opioid)
      set.seed(seed * 2)
      fx <- sample(0:1, size = n.noud, prob = c(1 - ini.NOUD.fx, ini.NOUD.fx), replace = T)
      pop.list[[t]]$fx[init_ppl$curr.state == "NODU"] <- fx
      m.tp <- trans.prob(pop.list[[t]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn <- init.Nx + n.nlx.yr / 12
    } else {
      pop.list[[t]] <- pop.list[[t - 1]]
      if (t %% 12 == 0) {
        OUD.fx <- min(ini.OUD.fx * (1 + gw.fx * min(floor((t - 1) / 12) + 1, 3)), 0.9)
        # determine fentanyl use among population who use opioids
        set.seed(seed)
        fx <- sample(0:1, size = n.opioid, prob = c(1 - OUD.fx, OUD.fx), replace = T)
        pop.list[[t]]$fx[init_ppl$curr.state != "NODU"] <- fx
        # # determine fentanyl exposure among population who use stimulants (non-opioid)
        # set.seed(seed*2)
        # fx         <- sample(0:1, size = n.noud, prob = c(1-ini.NOUD.fx, ini.NOUD.fx), replace = T)
        # pop.list[[t]]$fx[init_ppl$curr.state == "NODU"] <- fx
      }
      m.tp <- trans.prob(pop.list[[t - 1]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn <- n.nlx.mn * (1 - r.LossExp) + n.nlx.yr / 12
    }
    pop.list[[t]]$curr.state <- as.vector(samplev(probs = m.tp, m = 1)) # sample the next health state and store that state in matrix m.M
    ind.oustate.chg <- filter(pop.list[[t]], curr.state %in% v.oustate & OU.state != curr.state)$ind
    pop.list[[t]]$OU.state[ind.oustate.chg] <- pop.list[[t]]$curr.state[ind.oustate.chg]

    od_ppl <- pop.list[[t]][pop.list[[t]]$curr.state == "od", ]
    v.od[t] <- nrow(od_ppl)
    ou.pop.resid <- pop.list[[t]] %>% count(residence)

    decntree.out <- decisiotimesteptree(od_ppl, n.nlx = n.nlx.mn, ou.pop.resid, params, seed = seed + t)

    v.oddeath[t] <- sum(decntree.out[, "od.death"])
    v.oddeath.w[t] <- sum(decntree.out[decntree.out[, "wtns"] == 1, "od.death"])
    v.odpriv[t] <- sum(decntree.out[, "locpriv"])
    v.odpubl[t] <- v.od[t] - v.odpriv[t]
    v.deathpriv[t] <- sum(decntree.out[decntree.out[, "od.death"] == 1, "locpriv"])
    v.deathpubl[t] <- sum(decntree.out[, "od.death"] == 1) - v.deathpriv[t]
    v.nlxused[t] <- sum(decntree.out[, "nlx.used"])
    n.EMS <- sum(decntree.out[, "EMS"])
    n.hospcare <- sum(decntree.out[, "hospcare"])
    od_ppl$curr.state[decntree.out[, "od.death"] == 1] <- "dead"
    od_ppl$ever.od[decntree.out[, "od.death"] != 1] <- 1
    od_ppl$curr.state[decntree.out[, "inact"] == 1] <- "inact"
    od_ppl$curr.state[od_ppl$curr.state == "od"] <- od_ppl$OU.state[od_ppl$curr.state == "od"]

    m.oddeath.fx[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$fx == 1, ])
    m.oddeath.op[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU", ])
    m.oddeath.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU" & od_ppl$OU.state != "preb", ])
    m.oddeath.st[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "NODU", ])
    m.EDvisits[t] <- n.hospcare

    od.death.sum <- od_ppl[od_ppl$curr.state == "dead", ] %>% count(residence)
    for (dd in 1:nrow(od.death.sum)) {
      m.oddeath[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
    pop.list[[t]][od_ppl$ind, ] <- od_ppl
    cost.matrix[t, ] <- Costs(state = pop.list[[t]]$curr.state, OU.state = pop.list[[t]]$OU.state, nlx = sum(n.nlx.yr) / 12, count = list(n.EMS = n.EMS, n.hospcare = n.hospcare), params)

    pop.list[[t]]$age[pop.list[[t]]$curr.state != "dead"] <- pop.list[[t]]$init.age[pop.list[[t]]$curr.state != "dead"] + floor(t / 12) # update age for individuals that are still alive

    ## replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0)
    pop.list[[t]]$age[pop.list[[t]]$curr.state == "dead"] <- pop.list[[t]]$init.age[pop.list[[t]]$curr.state == "dead"]
    pop.list[[t]]$ever.od[pop.list[[t]]$curr.state == "dead"] <- 0
    pop.list[[t]]$OU.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state != "inact"] <- pop.list[[t]]$init.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state != "inact"]
    pop.list[[t]]$OU.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state == "inact"] <- pop.list[[t]]$OU.state[pop.list[[t]]$curr.state == "dead" & pop.list[[t]]$init.state == "inact"]
    pop.list[[t]]$curr.state[pop.list[[t]]$curr.state == "dead"] <- pop.list[[t]]$init.state[pop.list[[t]]$curr.state == "dead"]

    # cat('\r', paste(round(t/timesteps * 100), "% done", sep = " "))       # display the progress of the simulation
  } # end the loop for the time steps


  total.cost <- sum(cost.matrix[, "TotalCost"] * v.dwc) # total (discounted) cost

  if (PT.out == TRUE) {
    pop.trace <- pop.list
  } else {
    pop.trace <- NULL
  }

  results <- list(
    v.oddeath = v.oddeath, m.oddeath = m.oddeath, v.od = v.od,
    cost.matrix = cost.matrix, total.cost = total.cost, pop.trace = pop.trace, n.nlx.OEND.str = (n.nlx.mx.str - NxPharm.array[dim(NxPharm.array)[1], , ]), n.nlx.all.str = n.nlx.mx.str,
    m.oddeath.fx = m.oddeath.fx, m.oddeath.op = m.oddeath.op, m.oddeath.st = m.oddeath.st, m.oddeath.hr = m.oddeath.hr, m.EDvisits = m.EDvisits,
    v.odpriv = v.odpriv, v.odpubl = v.odpubl, v.deathpriv = v.deathpriv, v.deathpubl = v.deathpubl, v.nlxused = v.nlxused, v.oddeath.w = v.oddeath.w
  ) # store the results from the simulation in a list
  return(results) # return the results
} # end of the MicroSim function
