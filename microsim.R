#' Microsimulation for the PROFOUND model
#'
#' @description
#' `MicroSim()` runs the main model for naloxone distribution. It relies on the
#' decision tree function to simulate each individual's change of state and saves
#' population data including change of drug use status, non-fatal overdose, fatal
#' overdose, and cost.
#' 
#' `step()` handles one timestep of the model
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
#' `MicroSim()` returns a data frame of outcomes of the simulation
#' `step()` returns the results of one time step and the new population data frame
#'

source("transition_probability.R")
source("decision_tree.R")
source("cost_effectiveness.R")

MicroSim <- function(init_ppl, params, data, output, discount_rate, scenario = "SQ", seed = 1) {
  list2env(params, environment())
  # Find number of opioid and non-opioid users
  init_ppl.residence <- (init_ppl %>% count(residence))$n
  output <- data.frame(t = params$timesteps, scenario = scenario, v.od = rep(0, times = params$timesteps))

  # Pharmacy data
  NxPharm.mx <- data$nlx_data_pharm$pe[data$nlx_data_pharm$year >= (year_start - 1)] %*% t(init_ppl.residence / sum(init_ppl.residence))
  NxPharm.array <- array(0, dim = c(dim(NxPharm.mx)[1], 2, dim(NxPharm.mx)[2]))

  for (cc in 1:dim(NxPharm.mx)[1]) {
    NxPharm.array[cc, , ] <- round(rep(NxPharm.mx[cc, ], each = 2) * data$od_loc, 0)
  }
  nlx_array <- data$nx_oend[dimnames(data$nx_oend)[[1]] >= year_start, , ] + NxPharm.array[-1, , ]
  initial_nx <- data$nx_oend[dimnames(data$nx_oend)[[1]] == year_start - 1, , ] + NxPharm.array[1, , ]

  nlx_month <- nlx_array[dim(nlx_array)[1], , ]

  # Change naloxone based on scenario
  scenario <- scenarios[[scenario]]
  output$expansion <- scenario$expansion$val
  if (scenario$program$val) {
    avail_nlx <- nlx_month + evaluate_program() # TODO make this function work
  } else if (scenario$expansion$val > 1 && scenario$program$val) {
    avail_nlx <- data$nx_oend[dim(data$nx_oend)[1], , ] * scenario$expansion$val + NxPharm.array[dim(NxPharm.array)[1], , ]
  } else if (scenario$expansion$val > 1) {
    avail_nlx <- nlx_month + scenario$expansion$val
  } else {
    avail_nlx <- nlx_month
  }

  nlx_array <- abind(nlx_array, avail_nlx, along = 1)

  cost_discount <- rep(1 / (1 + discount_rate)^(0:(num_years - 1)), each = 12) # calculate the cost discount weight based on the discount rate

  # Create the population list to capture the state/attributes/costs for all individuals at each time point
  ppl_list <- list()
  ppl_list[[1]] <- init_ppl
  set.seed(seed) # set the seed for every individual for the random number generator
  # model step TODO: refactor into separate function
  for (t in 1:params$timesteps) {
    tmp <- step(t, output, nlx_array, ppl_list, data, seed, params, initial_nx)
    output <- tmp$output
    ppl_list <- tmp$ppl_list
  }

  total.cost <- sum(output$cost.matrix[, "Totalcost"] * cost_discount) # total (discounted) cost

  print("Saving results")
  return(output) # return the results
}

step <- function(t, output, nlx_array, ppl_list, data, seed, params, initial_nx) {
  # needs to update the ppl list, outputs
  output$t[t] <- t
  nx_avail_yr <- nlx_array[floor((t - 1) / 12) + 1, , ]
  if (t == 1) {
    set.seed(seed)
    num_opioid <- sum(ppl_list[[t]]$curr.state != "NODU")
    n.noud <- sum(ppl_list[[t]]$curr.state == "NODU")
    OUD.fx <- data$init_oud_fx
    fx <- sample(0:1, size = num_opioid, prob = c(1 - data$init_oud_fx, data$init_oud_fx), replace = T)
    ppl_list[[t]]$fx[ppl_list[[t]]$curr.state != "NODU"] <- fx
    # determine fentanyl use among population who use stimulants (non-opioid)
    set.seed(seed * 2)
    fx <- sample(0:1, size = n.noud, prob = c(1 - data$init_noud_fx, data$init_noud_fx), replace = T)
    ppl_list[[t]]$fx[ppl_list[[t]]$curr.state == "NODU"] <- fx
    trans_prob <- trans.prob(ppl_list[[t]], params, data) # calculate the transition probabilities at cycle t
    n.nlx.mn <- initial_nx + nx_avail_yr / 12
  } else {
    ppl_list[[t]] <- ppl_list[[t - 1]]
    if (t %% 12 == 0) {
      num_opioid <- sum(ppl_list[[t]]$curr.state != "NODU")
      n.noud <- sum(ppl_list[[t]]$curr.state == "NODU")
      OUD.fx <- min(data$init_oud_fx * (1 + data$fx_growth * min(floor((t - 1) / 12) + 1, 3)), 0.9)
      # determine fentanyl use among population who use opioids
      set.seed(seed)
      fx <- sample(0:1, size = num_opioid, prob = c(1 - OUD.fx, OUD.fx), replace = T)
      ppl_list[[t]]$fx[ppl_list[[t]]$curr.state != "NODU"] <- fx
      # determine fentanyl exposure among population who use stimulants (non-opioid)
      set.seed(seed * 2)
      fx <- sample(0:1, size = n.noud, prob = c(1 - data$init_noud_fx, data$init_noud_fx), replace = T)
      ppl_list[[t]]$fx[ppl_list[[t]]$curr.state == "NODU"] <- fx
    }
    trans_prob <- trans.prob(ppl_list[[t - 1]], params, data) # calculate the transition probabilities at cycle t
    n.nlx.mn <- initial_nx + nx_avail_yr / 12
    n.nlx.mn <- n.nlx.mn * (1 - data$r_loss_exp) + nx_avail_yr / 12
  }
  
  samples <- samplev(trans_prob, 1)

  # update health state
  ppl_list[[t]]$curr.state <- as.vector(samples)
  ind.oustate.chg <- filter(ppl_list[[t]], curr.state %in% params$oustate & OU.state != curr.state)$ind
  ppl_list[[t]]$OU.state[ind.oustate.chg] <- ppl_list[[t]]$curr.state[ind.oustate.chg]
  od_ppl <- ppl_list[[t]][ppl_list[[t]]$curr.state == "od", ]

  output$v.od <- nrow(od_ppl)
  ou.pop.resid <- ppl_list[[t]] %>% count(residence) %>% data.frame()
  rownames(ou.pop.resid) <- ou.pop.resid$residence
  decntree.out <- decision_tree(od_ppl, n.nlx.mn, ou.pop.resid, params, seed + t, data)

  output$oddeath[t] <- sum(decntree.out[, "od.death"])
  output$priv_od[t] <- sum(decntree.out[, "locpriv"])
  output$vpub_od[t] <- output$v.od[t] - output$priv_od[t]
  output$priv_death[t] <- sum(decntree.out[decntree.out[, "od.death"] == 1, "locpriv"])
  output$pub_death[t] <- sum(decntree.out[, "od.death"] == 1) - output$priv_death[t]
  output$nlxused[t] <- sum(decntree.out[, "nlx.used"])
  output$n.EMS[t] <- sum(decntree.out[, "EMS"])
  output$EDvisits[t] <- sum(decntree.out[, "hospcare"])
  od_ppl$curr.state[decntree.out[, "od.death"] == 1] <- "dead"
  od_ppl$ever_od[decntree.out[, "od.death"] != 1] <- 1
  od_ppl$curr.state[decntree.out[, "inact"] == 1] <- "inact"
  od_ppl$curr.state[od_ppl$curr.state == "od"] <- od_ppl$OU.state[od_ppl$curr.state == "od"]
  output$fx_deaths[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$fx == 1, ])
  output$oud_deaths[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU", ])
  output$hr_deaths[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU" & od_ppl$OU.state != "rx", ])
  output$stim_deaths[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "NODU", ])

  od.death.sum <- od_ppl[od_ppl$curr.state == "dead", ] %>% count(residence)
  if (nrow(od.death.sum) > 0) {
    for (dd in 1:nrow(od.death.sum)) {
      output[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
  }

  ppl_list[[t]][od_ppl$ind, ] <- od_ppl
  cost <- costs(state = ppl_list[[t]]$curr.state, OU.state = ppl_list[[t]]$OU.state, nlx = sum(nx_avail_yr) / 12, count = list(n.EMS = output$n.EMS, n.hospcare = output$n.hospcare), data)
  # should probably check if it's a new year and then increment it
  ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state != "dead"] + floor(t / 12) # update age for individuals that are still alive

  ## replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0)
  ppl_list[[t]]$age[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state == "dead"]
  ppl_list[[t]]$ever_od[ppl_list[[t]]$curr.state == "dead"] <- 0
  ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"]
  ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"]
  ppl_list[[t]]$curr.state[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead"]
  return(list(output = output, ppl_list = ppl_list))
}

