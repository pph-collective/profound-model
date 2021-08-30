#' Microsimulation for the PROFOUND model
#' 
#' @description 
#' `MicroSim()` runs the main model for naloxone distribution. It relies on the 
#' decision tree function to simulate each individual's change of state and saves 
#' population data including change of drug use status, non-fatal overdose, fatal
#' overdose, and cost.
#' 
#' @param init_ppl The initial population for the simulation.
#' @param params Model parameters.
#' @param data Empirical data to inform the model.
#' @param output The dataframe for model output.
#' @param discount_rate The rate of discounting for costs.
#' @param scenario The scenario to be run.
#' @param seed Random seed for simulation
#' 
#' @returns
#' outcomes of the simulation
#' 

###############################################################################################
#########################           Microsimulation        ####################################
###############################################################################################

###############################################################################################
####    Microsimulation to determine health states and number of overdoses                 ####
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


MicroSim <- function(init_ppl, params, data, output, discount_rate, scenario = "SQ", seed = 1) {
  # Arguments:
  # init_ppl:       matrix of initial states for individuals
  # params:    model parameters
  # timesteps:            total number of cycles to run the model
  # agent_states:        vector of health state names
  # discount.rate:  discount rate for costs
  # PT.out:         should the output include a Microsimulation trace? (default is TRUE)
  # scenario:            simulating strategy
  # seed:           starting seed number for random number generator (default is 1)
  # Makes use of:
  # trans.prob:     function for the estimation of transition probabilities
  # Costs:          function for the estimation of cost state values
  # decision_tree:  function for the decision tree module
  # TODO: actual docstring description
  list2env(params, environment())
  # Find number of opioid and non-opioid users
  num_opioid <- sum(init_ppl$curr.state != "NODU")
  n.noud <- sum(init_ppl$curr.state == "NODU")
  init_ppl.residence <- (init_ppl %>% count(residence))$n
  output <- data.frame(t = params$timesteps, scenario = scenario, v.od = rep(0, times = params$timesteps))

  # REVIEWED NxPharm is all data from pharmacy naloxone; only have overall number, so limited info
  NxPharm.mx <- data$NxDataPharm$pe[data$NxDataPharm$year >= (year_start - 1)] %*% t(init_ppl.residence / sum(init_ppl.residence))
  NxPharm.array <- array(0, dim = c(dim(NxPharm.mx)[1], 2, dim(NxPharm.mx)[2]))

  for (cc in 1:dim(NxPharm.mx)[1]) {
    NxPharm.array[cc, , ] <- round(rep(NxPharm.mx[cc, ], each = 2) * data$OD_loc, 0)
  }
  array.Nx <- data$NxOEND.array[dimnames(data$NxOEND.array)[[1]] >= year_start, , ] + NxPharm.array[-1, , ]
  initial_nx <- data$NxOEND.array[dimnames(data$NxOEND.array)[[1]] == year_start - 1, , ] + NxPharm.array[1, , ]

  n.nlx.mx.lst <- array.Nx[dim(array.Nx)[1], , ]

  scenario <- scenarios[[scenario]]
  output$expansion <- scenario$expansion$val
  if (scenario$program$val) {
    avail_nlx <- n.nlx.mx.lst + evaluate_program()  # TODO make this function work
  } else if (scenario$expansion$val > 1 && scenario$program$val) {
     avail_nlx <- data$NxOEND.array[dim(data$NxOEND.array)[1], , ] * scenario$expansion$val + NxPharm.array[dim(NxPharm.array)[1], , ]
  } else if (scenario$expansion$val > 1) {
    avail_nlx <- n.nlx.mx.lst + scenario$expansion$val
  }
  else {
    avail_nlx <- n.nlx.mx.lst 
  }

  array.Nx <- array.Nx <- abind(array.Nx, avail_nlx, along = 1)

  v.dwc <- rep(1 / (1 + discount_rate)^(0:(num_years - 1)), each = 12) # calculate the cost discount weight based on the discount rate

  # Create the population list to capture the state/attributes/costs for all individuals at each time point
  ppl_list <- list()
  set.seed(seed) # set the seed for every individual for the random number generator
  for (t in 1:params$timesteps) {
    output$t[t] <- t
    nx_avail_yr <- array.Nx[floor((t - 1) / 12) + 1, , ]
    if (t == 1) {
      ppl_list[[t]] <- init_ppl
      OUD.fx <- data$init_oud_fx
      # determine fentanyl use among population who use opioids
      set.seed(seed)

      fx <- sample(0:1, size = num_opioid, prob = c(1 - data$init_oud_fx, data$init_oud_fx), replace = T)
      ppl_list[[t]]$fx[init_ppl$curr.state != "NODU"] <- fx
      # determine fentanyl use among population who use stimulants (non-opioid)
      set.seed(seed * 2)
      fx <- sample(0:1, size = n.noud, prob = c(1 - data$ini.NOUD.fx, data$ini.NOUD.fx), replace = T)
      ppl_list[[t]]$fx[init_ppl$curr.state == "NODU"] <- fx
      m.tp <- trans.prob(ppl_list[[t]], params, data) # calculate the transition probabilities at cycle t
      n.nlx.mn <- initial_nx + nx_avail_yr / 12
    } else {
      ppl_list[[t]] <- ppl_list[[t - 1]]
      if (t %% 12 == 0) {
        OUD.fx <- min(data$init_oud_fx * (1 + data$gw.fx * min(floor((t - 1) / 12) + 1, 3)), 0.9)
        # determine fentanyl use among population who use opioids
        set.seed(seed)
        fx <- sample(0:1, size = num_opioid, prob = c(1 - OUD.fx, OUD.fx), replace = T)
        ppl_list[[t]]$fx[init_ppl$curr.state != "NODU"] <- fx
        # determine fentanyl exposure among population who use stimulants (non-opioid)
        set.seed(seed*2)
        fx         <- sample(0:1, size = n.noud, prob = c(1-data$ini.NOUD.fx, data$ini.NOUD.fx), replace = T)
        ppl_list[[t]]$fx[init_ppl$curr.state == "NODU"] <- fx
      }
      m.tp <- trans.prob(ppl_list[[t - 1]], params, data) # calculate the transition probabilities at cycle t
      n.nlx.mn <- n.nlx.mn * (1 - data$r.LossExp) + nx_avail_yr / 12
    }

    samples <- samplev(m.tp, 1)
    ppl_list[[t]]$curr.state <- as.vector(samples) # sample the next health state and store that state in matrix m.M
    ind.oustate.chg <- filter(ppl_list[[t]], curr.state %in% v.oustate & OU.state != curr.state)$ind
    ppl_list[[t]]$OU.state[ind.oustate.chg] <- ppl_list[[t]]$curr.state[ind.oustate.chg]

    od_ppl <- ppl_list[[t]][ppl_list[[t]]$curr.state == "od", ]

    output$v.od[t] <- nrow(od_ppl)

    ou.pop.resid <- ppl_list[[t]] %>% count(residence)
    decntree.out <- decision_tree(od_ppl, n.nlx = n.nlx.mn, ou.pop.resid, params, seed = seed + t)

    output$v.oddeath[t] <- sum(decntree.out[, "od.death"])
    output$v.odpriv[t] <- sum(decntree.out[, "locpriv"])
    output$v.odpubl[t] <- output$v.od[t] - output$v.odpriv[t]
    output$v.deathpriv[t] <- sum(decntree.out[decntree.out[, "od.death"] == 1, "locpriv"])
    output$v.deathpubl[t] <- sum(decntree.out[, "od.death"] == 1) - output$v.deathpriv[t]
    output$v.nlxused[t] <- sum(decntree.out[, "nlx.used"])
    output$n.EMS[t] <- sum(decntree.out[, "EMS"])
    output$EDvisits[t] <- sum(decntree.out[, "hospcare"])
    od_ppl$curr.state[decntree.out[, "od.death"] == 1] <- "dead"
    od_ppl$ever.od[decntree.out[, "od.death"] != 1] <- 1
    od_ppl$curr.state[decntree.out[, "inact"] == 1] <- "inact"
    od_ppl$curr.state[od_ppl$curr.state == "od"] <- od_ppl$OU.state[od_ppl$curr.state == "od"]

    output$m.oddeath.fx[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$fx == 1, ])
    output$m.oddeath.op[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU", ])
    output$m.oddeath.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU" & od_ppl$OU.state != "preb", ])
    output$m.oddeath.st[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "NODU", ])
    # output$m.EDvisits[t] <- output$n.hospcare

    od.death.sum <- od_ppl[od_ppl$curr.state == "dead", ] %>% count(residence)

    for (dd in 1:nrow(od.death.sum)) {
      output[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
      # output$m.oddeath[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }

    ppl_list[[t]][od_ppl$ind, ] <- od_ppl
    cost <- Costs(state = ppl_list[[t]]$curr.state, OU.state = ppl_list[[t]]$OU.state, nlx = sum(nx_avail_yr) / 12, count = list(n.EMS = output$n.EMS, n.hospcare = output$n.hospcare), data)
    ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state != "dead"] + floor(t / 12) # update age for individuals that are still alive
    
    ## replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0)
    ppl_list[[t]]$age[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state == "dead"]
    ppl_list[[t]]$ever.od[ppl_list[[t]]$curr.state == "dead"] <- 0
    ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"]
    ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"]
    ppl_list[[t]]$curr.state[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead"]
  } # end the loop for the time steps

  total.cost <- sum(output$cost.matrix[, "TotalCost"] * v.dwc) # total (discounted) cost

  print("Saving results")
  return(output) # return the results
} # end of the MicroSim function
