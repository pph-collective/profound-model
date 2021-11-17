#' Microsimulation for the PROFOUND model
#'
#' @description
#' `MicroSim()` runs the main model for naloxone distribution. It relies on the
#' decision tree function to simulate each individual's change of state and saves
#' population data including change of drug use status, non-fatal overdose, fatal
#' overdose, and cost.
#'
#' `step()` handles one time step of the model
#'
#' @param init_ppl The initial population for the simulation.
#' @param params Model parameters.
#' @param output The data frame for model output.
#' @param scenario The scenario to be run.
#'
#' @returns
#' `microsim()` returns a data frame of outcomes of the simulation
#' `step()` returns the results of one time step and the new population data frame
#'

source("transition_probability.R")
source("decision_tree.R")
source("cost_effectiveness.R")

microsim <- function(init_ppl, params, output, scenarios) {
  # Find number of agents by opioid use type
  num_opioid_nonrx <- sum(
    init_ppl$curr.state == "il.lr" | init_ppl$curr.state == "il.hr"
  ) + sum(init_ppl$curr.state == "relap" & init_ppl$OU.state != "preb")

  num_opioid_preb <- sum(init_ppl$curr.state == "preb") +
    sum(init_ppl$curr.state == "relap" & init_ppl$OU.state == "preb")
  # Count of population per residence
  init_ppl.residence <- (init_ppl %>% count(residence))$n

  NxPharm.mx <- NxDataPharm$pe[NxDataPharm$year >= (yr_start - 1)] %*%
    t(init_ppl.residence / sum(init_ppl.residence))

  NxPharm.array <- array(0, dim = c(dim(NxPharm.mx)[1], 2, dim(NxPharm.mx)[2]))

  for (cc in 1:dim(NxPharm.mx)[1]) {
    NxPharm.array[cc, , ] <- round(
      matrix(rep(NxPharm.mx[cc, ], each = 2),
        nrow = 2
      ) * c(1 - OD_loc_pub, OD_loc_pub), 0
    )
  }

  array.Nx <- params$oend[dimnames(oend)[[1]] >= yr_start, , ] +
    NxPharm.array[-1, , ]
  initial_nx <- params$oend[dimnames(oend)[[1]] == yr_start - 1, , ] +
    NxPharm.array[1, , ]
  n.nlx.mx.lst <- array.Nx[dim(array.Nx)[1], , ]
  
  if (strategy == "SQ") {
    n.nlx.mx.str <- n.nlx.mx.lst
  } else if (strategy == "expand") {
    # TODO what is exp.lv
    n.nlx.mx.str <- oend[dim(oend)[1], , ] * exp.lv + NxPharm.array[dim(NxPharm.array)[1], , ]
  } else if (strategy == "program") {
    # TODO what is pg.add
    n.nlx.mx.str <- n.nlx.mx.lst + pg.add
  }

  for (aa in 1:(num_years - dim(array.Nx)[1])) {
    array.Nx <- abind(array.Nx, n.nlx.mx.str, along = 1)
  }

  # cost discount weight
  v.dwc <- rep(1 / (1 + params$discount)^(0:(num_years - 1)), each = 12)

  # Create the population list to capture the state/attributes/costs for all individuals at each time point
  ppl_list <- list()
  set.seed(seed) # set the seed for every individual for the random number generator

  for (t in 1:timesteps) {
    nx_avail_yr <- array.Nx[floor((t - 1) / 12) + 1, , ]
    if (t == 1) {
      ppl_list[[t]] <- init_ppl
      OUD.fx <- ini.oud.fx
      params$p.preb2inact <- p.preb2inact.ini
      params$p.il.lr2inact <- p.il.lr2inact.ini
      params$p.il.hr2inact <- p.il.hr2inact.ini
      # determine fentanyl use among population who use non-prescription opioids (heroin)
      set.seed(seed)
      fx_nonpreb <- sample(0:1, size = num_opioid_nonpreb, prob = c(1 - OUD.fx, OUD.fx), replace = T)
      ppl_list[[t]]$fx[with(ppl_list[[t]], ind[curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")])] <- fx_nonpreb
      # determine fentanyl use among population who use prescription opioids
      set.seed(seed * 2)
      OUD.preb.fx <- OUD.fx * out.prebopioid # prevalence of fentanyl exposure is diluted by the portion outsourced opioids not from prescription
      fx_preb <- sample(0:1, size = num_opioid_preb, prob = c(1 - OUD.preb.fx, OUD.preb.fx), replace = T)
      ppl_list[[t]]$fx[with(ppl_list[[t]], ind[curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")])] <- fx_preb
      # determine fentanyl use among population who use stimulants (non-opioid)
      set.seed(seed * 3)
      fx_noud <- sample(0:1, size = num_noud, prob = c(1 - ini.NOUD.fx, ini.NOUD.fx), replace = T)
      ppl_list[[t]]$fx[ppl_list[[t]]$curr.state == "NODU"] <- fx_noud
      m.tp <- trans.prob(ppl_list[[t]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn <- initial_nx + nx_avail_yr / 12
    } else {
      ppl_list[[t]] <- ppl_list[[t - 1]]

      gw.2inact <- ifelse(t < (2019 - 2016 + 1) * 12, (1 + gw.m.2inact)^t, (1 + gw.m.2inact)^((2019 - 2016 + 1) * 12))
      params$p.preb2inact <- p.preb2inact.ini * gw.2inact
      params$p.il.lr2inact <- p.il.lr2inact.ini * gw.2inact
      params$p.il.hr2inact <- p.il.hr2inact.ini * gw.2inact

      if ((t - 1) %% 12 == 0) { # adjust and reassign fentanyl exposure every year
        OUD.fx <- min(ini.oud.fx * (1 + gw.fx * min(floor((t - 1) / 12) + 1, 3)), 0.9)
        num_opioid_nonpreb <- sum(ppl_list[[t]]$curr.state == "il.lr" | ppl_list[[t]]$curr.state == "il.hr") +
          sum(ppl_list[[t]]$curr.state == "relap" & ppl_list[[t]]$OU.state != "preb")
        num_opioid_preb <- sum(ppl_list[[t]]$curr.state == "preb") +
          sum(ppl_list[[t]]$curr.state == "relap" & ppl_list[[t]]$OU.state == "preb")
        num_noud <- sum(ppl_list[[t]]$curr.state == "NODU")

        # determine fentanyl use among population who use non-prescription opioids (heroin)
        set.seed(seed)
        fx_nonpreb <- sample(0:1, size = num_opioid_nonpreb, prob = c(1 - OUD.fx, OUD.fx), replace = T)
        ppl_list[[t]]$fx[with(ppl_list[[t]], ind[curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")])] <- fx_nonpreb
        # determine fentanyl use among population who use prescription opioids
        set.seed(seed * 2)
        fx_preb <- sample(0:1, size = num_opioid_preb, prob = c(1 - OUD.preb.fx, OUD.preb.fx), replace = T)
        ppl_list[[t]]$fx[with(ppl_list[[t]], ind[curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")])] <- fx_preb
        # # determine fentanyl exposure among population who use stimulants (non-opioid)
        # set.seed(seed*2)
        # fx         <- sample(0:1, size = n.noud, prob = c(1-ini.NOUD.fx, ini.NOUD.fx), replace = T)
        # ppl_list[[t]]$fx[init_ppl$curr.state == "NODU"] <- fx
      }
      m.tp <- trans.prob(ppl_list[[t - 1]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn <- n.nlx.mn * (1 - r.LossExp) + nx_avail_yr / 12
    }
    ppl_list[[t]]$curr.state <- as.vector(samplev(probs = m.tp, m = 1)) # sample the next health state and store that state in matrix m.M
    ind.oustate.chg <- filter(ppl_list[[t]], curr.state %in% v.oustate & OU.state != curr.state)$ind
    ppl_list[[t]]$OU.state[ind.oustate.chg] <- ppl_list[[t]]$curr.state[ind.oustate.chg]

    od_ppl <- ppl_list[[t]][ppl_list[[t]]$curr.state == "od", ]
    v.od[t] <- nrow(od_ppl)
    ou.pop.resid <- ppl_list[[t]] %>% count(residence)

    decntree.out <- decision_tree(od_ppl, n.nlx = n.nlx.mn, ou.pop.resid, params, seed = seed + t)

    v.oddeath[t] <- sum(decntree.out[, "od.death"])
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

    ## ADDED for od deaths stratified by population groups##
    m.oddeath.preb[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "preb", ])
    m.oddeath.il.lr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "il.lr", ])
    m.oddeath.il.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "il.hr", ])

    m.nlx.mn[t] <- sum(n.nlx.mn)
    ####

    od.death.sum <- od_ppl[od_ppl$curr.state == "dead", ] %>% count(residence)
    for (dd in 1:nrow(od.death.sum)) {
      m.oddeath[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
    ppl_list[[t]][od_ppl$ind, ] <- od_ppl
    cost.matrix[t, ] <- Costs(state = ppl_list[[t]]$curr.state, OU.state = ppl_list[[t]]$OU.state, nlx = sum(nx_avail_yr) / 12, count = list(n.EMS = n.EMS, n.hospcare = n.hospcare), params)

    ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state != "dead"] + floor(t / 12) # update age for individuals that are still alive

    ## replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0)
    ppl_list[[t]]$age[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state == "dead"]
    ppl_list[[t]]$ever.od[ppl_list[[t]]$curr.state == "dead"] <- 0
    ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"]
    ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"]
    ppl_list[[t]]$curr.state[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead"]

    # cat('\r', paste(round(t/timesteps * 100), "% done", sep = " "))       # display the progress of the simulation
  } # end the loop for the time steps


  total.cost <- sum(cost.matrix[, "TotalCost"] * v.dwc) # total (discounted) cost

  if (PT.out == TRUE) {
    pop.trace <- ppl_list
  } else {
    pop.trace <- NULL
  }
  print("Saving results")
  results <- list(
    v.oddeath = v.oddeath, m.oddeath = m.oddeath, v.od = v.od,
    cost.matrix = cost.matrix, total.cost = total.cost, pop.trace = pop.trace, n.nlx.OEND.str = (n.nlx.mx.str - NxPharm.array[dim(NxPharm.array)[1], , ]), avail_nlx = n.nlx.mx.str,
    m.oddeath.fx = m.oddeath.fx, m.oddeath.op = m.oddeath.op, m.oddeath.st = m.oddeath.st, m.oddeath.hr = m.oddeath.hr, m.EDvisits = m.EDvisits,
    v.odpriv = v.odpriv, v.odpubl = v.odpubl, v.deathpriv = v.deathpriv, v.deathpubl = v.deathpubl, v.nlxused = v.nlxused,
    m.oddeath.preb = m.oddeath.preb, m.oddeath.il.lr = m.oddeath.il.lr, m.oddeath.il.hr = m.oddeath.il.hr, m.nlx.mn = m.nlx.mn
  ) # store the results from the simulation in a list
  return(results) # return the results
} # end of the MicroSim function
