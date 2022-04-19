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

source("crosstab.R")
MicroSim <- function(init_ppl, params, timesteps, agent_states, discount.rate, PT.out = TRUE, strategy = "SQ", seed = 1) {
  # Arguments:
  # init_ppl:       matrix of initial states for individuals
  # params:         model parameters
  # timesteps:      total number of cycles to run the model
  # agent_states:   vector of health state names
  # discount.rate:  discount rate for costs
  # PT.out:         should the output include a Microsimulation trace? (default is TRUE)
  # strategy:       simulating strategy
  # seed:           starting seed number for random number generator (default is 1)
  # Makes use of:
  # trans.prob:     function for the estimation of transition probabilities
  # Costs:          function for the estimation of cost state values
  # decision_tree:  function for the decision tree module
  # TODO: actual docstring description
  list2env(params, environment())
  
  init_ppl.residence <- (init_ppl %>% count(residence))$n
  # REVIEWED NxPharm is all data from pharmacy naloxone; only have overall number, so limited info
  NxPharm.matrix <- round(NxDataPharm$pe[NxDataPharm$year >= (yr_start - 1)] %*% t(init_ppl.residence / sum(init_ppl.residence)), 0)
  rownames(NxPharm.matrix) <- NxDataPharm$year
  colnames(NxPharm.matrix) <- v.region

  matrix.Nx <- list(OEND = NxOEND.matrix[dimnames(NxOEND.matrix)[[1]] >= yr_start, ], Pharm = NxPharm.matrix[dimnames(NxPharm.matrix)[[1]] >= yr_start, ])
  initial.Nx <- list(OEND = NxOEND.matrix[dimnames(NxOEND.matrix)[[1]] == yr_start - 1, ], Pharm = NxPharm.matrix[dimnames(NxPharm.matrix)[[1]] == yr_start - 1, ])

  n.nlx.lst <- list(OEND = NxOEND.matrix[dim(NxOEND.matrix)[1], ], Pharm = NxPharm.matrix[dim(NxPharm.matrix)[1], ])
  if (strategy == "SQ") {
    n.nlx.str <- n.nlx.lst
  } else if (strategy == "expand") {
    n.nlx.str <- list(OEND = n.nlx.lst$OEND * exp.lv,
                      Pharm = n.nlx.lst$Pharm)
  } else if (strategy == "program") {
    n.nlx.str <- list(OEND = n.nlx.lst$OEND + pg.add,
                      Pharm = n.nlx.lst$Pharm)
    # n.nlx.str <- n.nlx.lst + pg.add
  } else if (grepl("10K", strategy, fixed = T)) {
    matrix.Nx$OEND <- cbind(matrix.Nx$OEND, tenk = numeric(dim(matrix.Nx$OEND)[1]))
    initial.Nx$OEND <- cbind(initial.Nx$OEND, tenk = 0)
    n.nlx.str <- list(OEND = cbind(n.nlx.lst$OEND, tenk = expand.kits),
                      Pharm = n.nlx.lst$Pharm)
  }

  for (aa in 1:(num_years - dim(matrix.Nx$OEND)[1])) {
    matrix.Nx$OEND <- rbind(matrix.Nx$OEND, n.nlx.str$OEND)
    matrix.Nx$Pharm <- rbind(matrix.Nx$Pharm, n.nlx.str$Pharm)
  }

  v.dwc <- rep(1 / (1 + discount.rate)^(0:(num_years - 1)), each = 12) # calculate the cost discount weight based on the discount rate

  # Create the population list to capture the state/attributes/costs for all individuals at each time point
  ppl_list <- list()
  
  for (t in 1:timesteps) {
    nx_avail_yr <- list(OEND = matrix.Nx$OEND[floor((t - 1) / 12) + 1, ],
                        Pharm = matrix.Nx$Pharm[floor((t - 1) / 12) + 1, ])
    if (t == 1) {
      ppl_list[[t]] <- init_ppl
      OUD.fx <- ini.oud.fx
      NOUD.fx <- ini.NOUD.fx
      params$p.preb2inact <- p.preb2inact.ini
      params$p.il.lr2inact <- p.il.lr2inact.ini
      params$p.il.hr2inact <- p.il.hr2inact.ini
      # determine number of opioid and non-opioid users who are susceptible to fentanyl
       #find the index of individuals using heroin (non-prescription opioids), all of which are at risk for fentanyl exposure, including those in the relapsed state
      ind_opioid_nonpreb <- with(init_ppl, ind[curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")])
       #find the index of individuals using prescription opioids, only a portion of them outsourced their opioids not from prescription (at risk for fentanyl exposure), including those in the relapsed state
      ind_opioid_preb <- with(init_ppl, ind[curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")])
       #find the index of individuals using stimulants (non-opioid drug use), all are at risk for fentanyl exposure
      ind_NODU <- with(init_ppl, ind[curr.state == "NODU"])
      # determine fentanyl exposure among population who use non-prescription opioids (heroin)
      set.seed(seed)
      ind_fx_nonpreb <- sample(ind_opioid_nonpreb, size = round(length(ind_opioid_nonpreb) * OUD.fx, 0))
      ppl_list[[t]]$fx[ind_fx_nonpreb] <- 1
      # determine fentanyl exposure among population who use prescription opioids
      OUD.preb.fx <- OUD.fx * out.prebopioid   #prevalence of fentanyl exposure is diluted by the portion outsourced opioids not from prescription
      set.seed(seed*2)
      ind_fx_preb <- sample(ind_opioid_preb, size = round(length(ind_opioid_preb) * OUD.preb.fx, 0))
      ppl_list[[t]]$fx[ind_fx_preb] <- 1
      # determine fentanyl exposure among population who use stimulants (non-opioid)
      set.seed(seed*3)
      ind_fx_NODU <- sample(ind_NODU, size = round(length(ind_NODU) * NOUD.fx, 0))
      ppl_list[[t]]$fx[ind_fx_NODU] <- 1
      m.tp <- trans.prob(ppl_list[[t]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn <- list(OEND = initial.Nx$OEND + nx_avail_yr$OEND / 12,
                       Pharm = initial.Nx$Pharm + nx_avail_yr$Pharm / 12)
    } else {
      ppl_list[[t]] <- ppl_list[[t - 1]]
      gw.2inact <- ifelse(t < (2019-2016 + 1)*12, (1+gw.m.2inact)^t, (1+gw.m.2inact)^((2019-2016 + 1)*12))
      OUD.fx <- min(ini.oud.fx * ifelse(t < (2019-2016 + 1)*12, (1+gw.fx)^t, (1+gw.fx)^((2019-2016 + 1)*12)), 0.9)
      NOUD.fx <- min(ini.NOUD.fx * ifelse(t < (2019-2016 + 1)*12, (1+gw.fx)^t, (1+gw.fx)^((2019-2016 + 1)*12)), 0.9)

      params$p.preb2inact <- p.preb2inact.ini * gw.2inact
      params$p.il.lr2inact <- p.il.lr2inact.ini * gw.2inact
      params$p.il.hr2inact <- p.il.hr2inact.ini * gw.2inact
      ppl_list[[t]]$fx <- 0
      # determine index of opioid and non-opioid users who are susceptible to fentanyl
      ind_opioid_nonpreb <- with(ppl_list[[t]], ind[curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")])
      ind_opioid_preb <- with(ppl_list[[t]], ind[curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")])
      ind_NODU <- with(ppl_list[[t]], ind[curr.state == "NODU"])
      # determine fentanyl use among population who use non-prescription opioids (heroin)
      set.seed(seed)
      ind_fx_nonpreb <- sample(ind_opioid_nonpreb, size = round(length(ind_opioid_nonpreb) * OUD.fx, 0))
      ppl_list[[t]]$fx[ind_fx_nonpreb] <- 1

      # determine fentanyl use among population who use prescription opioids
      OUD.preb.fx <- OUD.fx * out.prebopioid   #prevalence of fentanyl exposure is diluted by the portion outsourced opioids not from prescription
      set.seed(seed*2)
      ind_fx_preb <- sample(ind_opioid_preb, size = round(length(ind_opioid_preb) * OUD.preb.fx, 0))
      ppl_list[[t]]$fx[ind_fx_preb] <- 1

      # determine fentanyl exposure among population who use stimulants (non-opioid)
      set.seed(seed*3)
      ind_fx_NODU <- sample(ind_NODU, size = round(length(ind_NODU) * NOUD.fx, 0))
      ppl_list[[t]]$fx[ind_fx_NODU] <- 1
      m.tp <- trans.prob(ppl_list[[t]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn$OEND <- n.nlx.mn$OEND * (1 - r.LossExp) + nx_avail_yr$OEND / 12
      n.nlx.mn$Pharm <- n.nlx.mn$Pharm * (1 - r.LossExp) + nx_avail_yr$Pharm / 12
    }

    set.seed(seed) # set the seed for every individual for the random number generator
    ppl_list[[t]]$curr.state <- as.vector(samplev(probs = m.tp, m = 1)) # sample the next health state
    ind.oustate.chg <- filter(ppl_list[[t]], curr.state %in% v.oustate & OU.state != curr.state)$ind
    ppl_list[[t]]$OU.state[ind.oustate.chg] <- ppl_list[[t]]$curr.state[ind.oustate.chg]
    od_ppl <- ppl_list[[t]][ppl_list[[t]]$curr.state == "od", ]
    v.od[t] <- nrow(od_ppl)
    # ou.pop.resid <- ppl_list[[t]] %>% count(residence)
    crosstab_drug_resid <- crosstab(ppl_list[[t]], row.vars = "residence", col.vars = "OU.state", type = "f")$table
    crosstab_drug_resid <- crosstab_drug_resid[-nrow(crosstab_drug_resid), -ncol(crosstab_drug_resid)]
    # Naloxone availability algorithm to determine the provability of naloxone used in a witnessed overdose
    p.nlx.avail.mx <- nlx.avail.algm(n.nlx = n.nlx.mn, crosstab_drug_resid, OD_loc_pub, nlx.adj, cap, eff.pharNlx, strategy, t)
    decntree.out <- decision_tree(od_ppl, p.nlx.avail.mx, params, strategy, t, seed = seed + t)

    v.oddeath[t] <- sum(decntree.out[, "od.death"])
    v.oddeath.w[t] <- sum(decntree.out[decntree.out[ , "wtns"] == 1 , "od.death"])
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
    # m.oddeath.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU" & od_ppl$OU.state != "preb", ])
    m.oddeath.st[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "NODU", ])
    m.EDvisits[t] <- n.hospcare
    
    ##ADDED for od deaths stratified by population groups##
    m.oddeath.preb[t]  <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state =="preb", ])
    m.oddeath.il.lr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state =="il.lr", ])
    m.oddeath.il.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state =="il.hr", ])
    
    m.nlx.mn.OEND[t] <- sum(n.nlx.mn$OEND)
    m.nlx.mn.Pharm[t] <- sum(n.nlx.mn$Pharm)
    ####

    od.death.sum <- od_ppl[od_ppl$curr.state == "dead", ] %>% count(residence)
    for (dd in 1:nrow(od.death.sum)) {
      m.oddeath[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
    ppl_list[[t]][od_ppl$ind, ] <- od_ppl
    cost.matrix[t, ] <- Costs(state = ppl_list[[t]]$curr.state, OU.state = ppl_list[[t]]$OU.state, nlx = sum(unlist(nx_avail_yr)) / 12, count = list(n.EMS = n.EMS, n.hospcare = n.hospcare), params)

    ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] <- ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] + ifelse(t %% 12 == 0, 1, 0) # update age for individuals that are still alive

    ## replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0)
    ppl_list[[t]]$age[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state == "dead"]
    ppl_list[[t]]$ever.od[ppl_list[[t]]$curr.state == "dead"] <- 0
    ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"]
    ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"]
    ppl_list[[t]]$curr.state[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead"]
    # cat('\r', paste(round(t/timesteps * 100), "% done", sep = " "))       # display the progress of the simulation
  } # end the loop for the time steps


  total.cost <- sum(cost.matrix[, "TotalCost"] * v.dwc) # total (discounted) cost

  if (PT.out == TRUE) {
    pop.trace <- ppl_list
  } else {
    pop.trace <- NULL
  }
  if (grepl("10K", strategy)){
    n.nlx.OEND.str.base <- n.nlx.str$OEND[, -ncol(n.nlx.str$OEND)]
    n.nlx.OEND.str.add.total <- n.nlx.str$OEND[, ncol(n.nlx.str$OEND)]
    if (strategy == "SSP_10K"){
      SSP_drug_resid <- crosstab_drug_resid[, c("il.hr", "NODU")]*rep(c(1, 0.133), each = nrow(crosstab_drug_resid))
      n.nlx.OEND.str <- round(rowSums(n.nlx.OEND.str.add.total * SSP_drug_resid / sum(SSP_drug_resid)), 0) + n.nlx.OEND.str.base
    # } else if (strategy == "StreetOutreach_10K"){
    #   StreetOutreach_drug_resid <- crosstab_drug_resid[, c("il.lr", "il.hr", "NODU")]
    #   n.nlx.OEND.str <- round(rowSums(n.nlx.OEND.str.add.total * StreetOutreach_drug_resid / sum(StreetOutreach_drug_resid)), 0) + n.nlx.OEND.str.base
    } else if (strategy == "MailEvent_10K"){
      MailEvent_drug_resid <- crosstab_drug_resid
      n.nlx.OEND.str <- round(rowSums(n.nlx.OEND.str.add.total * MailEvent_drug_resid / sum(MailEvent_drug_resid)), 0) + n.nlx.OEND.str.base
    } else if (strategy == "Healthcare_10K"){
      Healthcare_drug_resid <- crosstab_drug_resid[, c("preb")]
      n.nlx.OEND.str <- round(n.nlx.OEND.str.add.total * Healthcare_drug_resid / sum(Healthcare_drug_resid), 0) + n.nlx.OEND.str.base
    }
  } else {
    n.nlx.OEND.str <- n.nlx.str$OEND
  }
  print("Saving results")
  results <- list(
    v.oddeath = v.oddeath, m.oddeath = m.oddeath, v.od = v.od,
    cost.matrix = cost.matrix, total.cost = total.cost, pop.trace = pop.trace,  n.nlx.OEND.str = n.nlx.OEND.str,
    m.oddeath.fx = m.oddeath.fx, m.oddeath.op = m.oddeath.op, m.oddeath.st = m.oddeath.st, m.oddeath.hr = m.oddeath.hr, m.EDvisits = m.EDvisits,
    v.odpriv = v.odpriv, v.odpubl = v.odpubl, v.deathpriv = v.deathpriv, v.deathpubl = v.deathpubl, v.nlxused = v.nlxused, v.oddeath.w = v.oddeath.w,
    m.oddeath.preb = m.oddeath.preb, m.oddeath.il.lr = m.oddeath.il.lr, m.oddeath.il.hr = m.oddeath.il.hr, m.nlx.mn.OEND = m.nlx.mn.OEND, m.nlx.mn.Pharm = m.nlx.mn.Pharm, 
    end.pop.state = table(ppl_list[[t]]$curr.state)
  ) # store the results from the simulation in a list
  return(results) # return the results
} # end of the MicroSim function
