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
  
  pop.summary <- init_ppl %>% count(residence, race)
  ## ADJUSTEMENT FOR CATCHMENTS WITH 0 POPULATION FOR ONE RACE
  for (pp in 1:length(unique(pop.summary$residence))){
    current.res <- unique(pop.summary$residence)[pp]
    current.pop.cat <-  subset(pop.summary, residence == current.res)
    if(nrow(current.pop.cat)<3){
      if(!"black" %in% current.pop.cat$race){
        pop.summary <- rbind(pop.summary, c(current.res, "black", 1))
      } else if (!"hisp" %in% current.pop.cat$race){
        pop.summary <- rbind(pop.summary, c(current.res, "hisp", 1))
      }
    }
  }
  pop.summary <- pop.summary[order(pop.summary$residence, pop.summary$race), ]
  pop.summary$n <- as.numeric(pop.summary$n)
  #### END ADJUSTMENT ####
  init_ppl.residence <- cbind(pop.summary$n[pop.summary$race=="white"], pop.summary$n[pop.summary$race=="black"], pop.summary$n[pop.summary$race=="hisp"])
  rownames(init_ppl.residence) <- v.region
  colnames(init_ppl.residence) <- c("white", "black", "hisp")
  # REVIEWED NxPharm is all data from pharmacy naloxone; only have overall number, so limited info
  NxPharm.array <- array(0, dim = c(length(NxDataPharm$year), length(v.region), 3))
  dimnames(NxPharm.array)[[1]] <- NxDataPharm$year
  dimnames(NxPharm.array)[[2]] <- v.region
  dimnames(NxPharm.array)[[3]] <- c("white", "black", "hisp")
  NxPharm.array[,,"white"] <- round(NxDataPharm$white[NxDataPharm$year >= (yr_start - 1)] %*% t(init_ppl.residence[,"white"] / sum(init_ppl.residence[,"white"])), 0)
  NxPharm.array[,,"black"] <- round(NxDataPharm$black[NxDataPharm$year >= (yr_start - 1)] %*% t(init_ppl.residence[,"black"] / sum(init_ppl.residence[,"black"])), 0)
  NxPharm.array[,,"hisp"]  <- round(NxDataPharm$hisp[NxDataPharm$year >= (yr_start - 1)] %*% t(init_ppl.residence[,"hisp"] / sum(init_ppl.residence[,"hisp"])), 0)
  NxOEND.array <- NxOEND.array[,,-4]
  array.Nx   <- list(OEND = NxOEND.array[dimnames(NxOEND.array)[[1]] >= yr_start, ,], Pharm = NxPharm.array[dimnames(NxPharm.array)[[1]] >= yr_start, ,])
  initial.Nx <- list(OEND = NxOEND.array[dimnames(NxOEND.array)[[1]] == yr_start - 1, ,], Pharm = NxPharm.array[dimnames(NxPharm.array)[[1]] == yr_start - 1, ,])
  n.nlx.lst  <- list(OEND = NxOEND.array[dim(NxOEND.array)[1], ,], Pharm = NxPharm.array[dim(NxPharm.array)[1], ,])
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
  } else if (grepl("ratio", strategy, fixed = T)) {
    if (grepl("NT", strategy, fixed = T)){
      current.ratio <- sum(n.nlx.lst$OEND) / sum(ODdeath2020[-4])
      n.nlx.str <- list(OEND = round(n.nlx.lst$OEND * (Nx2ODratio / current.ratio),0),
                        Pharm = n.nlx.lst$Pharm)
    } else if (grepl("EQ", strategy, fixed = T)){
      current.ratio <- colSums(n.nlx.lst$OEND) / ODdeath2020[-4]
      OEND.temp <- n.nlx.lst$OEND
      OEND.temp[,"white"] <- round(as.numeric(Nx2ODratio / current.ratio[1]) * OEND.temp[,"white"],0)
      OEND.temp[,"black"] <- round(as.numeric(Nx2ODratio / current.ratio[2]) * OEND.temp[,"black"],0)
      OEND.temp[,"hisp"]  <- round(as.numeric(Nx2ODratio / current.ratio[3]) * OEND.temp[,"hisp"],0)
      n.nlx.str <- list(OEND = OEND.temp,
                        Pharm = n.nlx.lst$Pharm)
    }
  }

  for (aa in 1:(num_years - dim(array.Nx$OEND)[1])) {
    array.Nx$OEND  <- abind(array.Nx$OEND, n.nlx.str$OEND, along =1)
    array.Nx$Pharm <- abind(array.Nx$Pharm, n.nlx.str$Pharm, along =1)
  }

  v.dwc <- rep(1 / (1 + discount.rate)^(1:3), each = 12) # calculate the cost and qaly discount weight based on the discount rate

  # Create the population list to capture the state/attributes/costs for all individuals at each time point
  ppl_list <- list()
  
  for (t in 1:timesteps) {
    nx_avail_yr <- list(OEND = array.Nx$OEND[floor((t - 1) / 12) + 1, ,],
                        Pharm = array.Nx$Pharm[floor((t - 1) / 12) + 1, ,])
    if (t == 1) {
      ppl_list[[t]] <- init_ppl
      OUD.fx <- ini.oud.fx
      NODU.fx <- ini.NODU.fx
      params$p.preb2inact <- p.preb2inact.ini
      params$p.il.lr2inact <- p.il.lr2inact.ini
      params$p.il.hr2inact <- p.il.hr2inact.ini
      # determine number of opioid and non-opioid users who are susceptible to fentanyl
       #find the index of individuals using heroin (non-prescription opioids), all of which are at risk for fentanyl exposure, including those in the relapsed state
      ind_opioid_nonpreb_white <- with(init_ppl, ind[(curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")) & race =="white"])
      ind_opioid_nonpreb_black <- with(init_ppl, ind[(curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")) & race =="black"])
      ind_opioid_nonpreb_hisp  <- with(init_ppl, ind[(curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")) & race =="hisp"])
       #find the index of individuals using prescription opioids, only a portion of them outsourced their opioids not from prescription (at risk for fentanyl exposure), including those in the relapsed state
      ind_opioid_preb_white <- with(init_ppl, ind[(curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")) & race == "white"])
      ind_opioid_preb_black <- with(init_ppl, ind[(curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")) & race == "black"])
      ind_opioid_preb_hisp  <- with(init_ppl, ind[(curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")) & race == "hisp"])
       #find the index of individuals using stimulants (non-opioid drug use), all are at risk for fentanyl exposure
      ind_NODU_white <- with(init_ppl, ind[curr.state == "NODU" & race == "white"])
      ind_NODU_black <- with(init_ppl, ind[curr.state == "NODU" & race == "black"])
      ind_NODU_hisp  <- with(init_ppl, ind[curr.state == "NODU" & race == "hisp"])
      # assign fentanyl exposure among population who use non-prescription opioids (heroin)
      set.seed(seed)
      ind_fx_nonpreb_white <- sample(ind_opioid_nonpreb_white, size = round(length(ind_opioid_nonpreb_white) * OUD.fx["white"], 0))
      ind_fx_nonpreb_black <- sample(ind_opioid_nonpreb_black, size = round(length(ind_opioid_nonpreb_black) * OUD.fx["black"], 0))
      ind_fx_nonpreb_hisp  <- sample(ind_opioid_nonpreb_hisp,  size = round(length(ind_opioid_nonpreb_hisp) * OUD.fx["hisp"], 0))
      ppl_list[[t]]$fx[which(ppl_list[[t]]$ind %in% c(ind_fx_nonpreb_white, ind_fx_nonpreb_black, ind_fx_nonpreb_hisp))] <- 1
      # assign fentanyl exposure among population who use prescription opioids
      OUD.preb.fx <- OUD.fx * out.prebopioid   #prevalence of fentanyl exposure is diluted by the portion outsourced opioids not from prescription
      set.seed(seed*2)
      ind_fx_preb_white <- sample(ind_opioid_preb_white, size = round(length(ind_opioid_preb_white) * OUD.preb.fx["white"], 0))
      ind_fx_preb_black <- sample(ind_opioid_preb_black, size = round(length(ind_opioid_preb_black) * OUD.preb.fx["black"], 0))
      ind_fx_preb_hisp  <- sample(ind_opioid_preb_hisp, size = round(length(ind_opioid_preb_hisp) * OUD.preb.fx["hisp"], 0))
      ppl_list[[t]]$fx[which(ppl_list[[t]]$ind %in% c(ind_fx_preb_white, ind_fx_preb_black, ind_fx_preb_hisp))] <- 1
      # assign fentanyl exposure among population who use stimulants (non-opioid)
      set.seed(seed*3)
      ind_fx_NODU_white <- sample(ind_NODU_white, size = round(length(ind_NODU_white) * NODU.fx["white"], 0))
      ind_fx_NODU_black <- sample(ind_NODU_black, size = round(length(ind_NODU_black) * NODU.fx["black"], 0))
      ind_fx_NODU_hisp  <- sample(ind_NODU_hisp, size = round(length(ind_NODU_hisp) * NODU.fx["hisp"], 0))
      ppl_list[[t]]$fx[which(ppl_list[[t]]$ind %in% c(ind_fx_NODU_white, ind_fx_NODU_black, ind_fx_NODU_hisp))] <- 1
      m.tp <- trans.prob(ppl_list[[t]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn <- list(OEND = initial.Nx$OEND + nx_avail_yr$OEND / 12,
                       Pharm = initial.Nx$Pharm + nx_avail_yr$Pharm / 12)
    } else {
      ppl_list[[t]] <- ppl_list[[t - 1]]
      
      if (t <= (2018-2016 + 1)*12){
        gw.2inact <- (1+gw.m.2inact)^t
        OUD.fx <- ini.oud.fx * (1+gw.fx)^t; if (any(OUD.fx > 0.95)){OUD.fx[which(OUD.fx > 0.95)] = 0.95}
        NODU.fx <- ini.NODU.fx * (1+gw.fx)^t; if (any(NODU.fx > 0.95)){NODU.fx[which(NODU.fx > 0.95)] = 0.95}
        # NODU.fx <- ini.NODU.fx + floor((t - 1) / 12)*gw.NODU.fx.ab.yr; if (any(NODU.fx > 0.9)){NODU.fx[which(NODU.fx > 0.9)] = 0.9}
      } else if (t <= (2019-2016 + 1)*12 + 2) {
        gw.2inact <- (1+gw.m.2inact)^t
        OUD.fx  <- ini.oud.fx * (1+gw.fx)^t; if (any(OUD.fx > 0.95)){OUD.fx[which(OUD.fx > 0.95)] = 0.95}
        NODU.fx  <- ini.NODU.fx * (1+gw.fx)^t; if (any(NODU.fx > 0.95)){NODU.fx[which(NODU.fx > 0.95)] = 0.95}
        # NODU.fx <- ini.NODU.fx + floor((t - 1) / 12)*gw.NODU.fx.ab.yr; if (any(NODU.fx > 0.95)){NODU.fx[which(NODU.fx > 0.95)] = 0.95}
      } else {
        gw.2inact <- (1+gw.m.2inact)^((2019-2016 + 1)*12+2)* (1- c(covid.rd.2inact.white, covid.rd.2inact.black, covid.rd.2inact.hisp))
        OUD.fx  <- ini.oud.fx * (1+gw.fx)^t; if (any(OUD.fx > 0.95)){OUD.fx[which(OUD.fx > 0.95)] = 0.95}
        NODU.fx  <- ini.NODU.fx * (1+gw.fx)^t; if (any(NODU.fx > 0.95)){NODU.fx[which(NODU.fx > 0.95)] = 0.95}
        # NODU.fx <- (ini.NODU.fx + (2019 - 2016) *gw.NODU.fx.ab.yr) * covid.NODU.fx + floor((t - 1) / 12 - 3)*gw.NODU.fx.ab.yr; if (any(NODU.fx > 0.9)){NODU.fx[which(NODU.fx > 0.9)] = 0.9}
      }
      params$p.preb2inact  <- p.preb2inact.ini * gw.2inact
      params$p.il.lr2inact <- p.il.lr2inact.ini * gw.2inact
      params$p.il.hr2inact <- p.il.hr2inact.ini * gw.2inact
      ppl_list[[t]]$fx <- 0
      # determine index of opioid and non-opioid users who are susceptible to fentanyl
      # ind_opioid_nonpreb <- with(ppl_list[[t]], ind[curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")])
      # ind_opioid_preb <- with(ppl_list[[t]], ind[curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")])
      # ind_NODU <- with(ppl_list[[t]], ind[curr.state == "NODU"])
      #find the index of individuals using heroin (non-prescription opioids), all of which are at risk for fentanyl exposure, including those in the relapsed state
      ind_opioid_nonpreb_white <- with(ppl_list[[t]], ind[(curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")) & race =="white"])
      ind_opioid_nonpreb_black <- with(ppl_list[[t]], ind[(curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")) & race =="black"])
      ind_opioid_nonpreb_hisp  <- with(ppl_list[[t]], ind[(curr.state == "il.lr" | curr.state == "il.hr" | (curr.state == "relap" & OU.state != "preb")) & race =="hisp"])
      #find the index of individuals using prescription opioids, only a portion of them outsourced their opioids not from prescription (at risk for fentanyl exposure), including those in the relapsed state
      ind_opioid_preb_white <- with(ppl_list[[t]], ind[(curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")) & race == "white"])
      ind_opioid_preb_black <- with(ppl_list[[t]], ind[(curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")) & race == "black"])
      ind_opioid_preb_hisp  <- with(ppl_list[[t]], ind[(curr.state == "preb" | (curr.state == "relap" & OU.state == "preb")) & race == "hisp"])
      #find the index of individuals using stimulants (non-opioid drug use), all are at risk for fentanyl exposure
      ind_NODU_white <- with(ppl_list[[t]], ind[curr.state == "NODU" & race == "white"])
      ind_NODU_black <- with(ppl_list[[t]], ind[curr.state == "NODU" & race == "black"])
      ind_NODU_hisp  <- with(ppl_list[[t]], ind[curr.state == "NODU" & race == "hisp"])
      
      # determine fentanyl use among population who use non-prescription opioids (heroin)
      # set.seed(seed)
      set.seed(seed + t)
      ind_fx_nonpreb_white <- sample(ind_opioid_nonpreb_white, size = round(length(ind_opioid_nonpreb_white) * OUD.fx["white"], 0))
      ind_fx_nonpreb_black <- sample(ind_opioid_nonpreb_black, size = round(length(ind_opioid_nonpreb_black) * OUD.fx["black"], 0))
      ind_fx_nonpreb_hisp  <- sample(ind_opioid_nonpreb_hisp,  size = round(length(ind_opioid_nonpreb_hisp) * OUD.fx["hisp"], 0))
      ppl_list[[t]]$fx[which(ppl_list[[t]]$ind %in% c(ind_fx_nonpreb_white, ind_fx_nonpreb_black, ind_fx_nonpreb_hisp))] <- 1

      # determine fentanyl use among population who use prescription opioids
      OUD.preb.fx <- OUD.fx * out.prebopioid   #prevalence of fentanyl exposure is diluted by the portion outsourced opioids not from prescription
      set.seed((seed + t)*2)
      ind_fx_preb_white <- sample(ind_opioid_preb_white, size = round(length(ind_opioid_preb_white) * OUD.preb.fx["white"], 0))
      ind_fx_preb_black <- sample(ind_opioid_preb_black, size = round(length(ind_opioid_preb_black) * OUD.preb.fx["black"], 0))
      ind_fx_preb_hisp  <- sample(ind_opioid_preb_hisp, size = round(length(ind_opioid_preb_hisp) * OUD.preb.fx["hisp"], 0))
      ppl_list[[t]]$fx[which(ppl_list[[t]]$ind %in% c(ind_fx_preb_white, ind_fx_preb_black, ind_fx_preb_hisp))] <- 1

      # determine fentanyl exposure among population who use stimulants (non-opioid)
      set.seed((seed + t)*3)
      ind_fx_NODU_white <- sample(ind_NODU_white, size = round(length(ind_NODU_white) * NODU.fx["white"], 0))
      ind_fx_NODU_black <- sample(ind_NODU_black, size = round(length(ind_NODU_black) * NODU.fx["black"], 0))
      ind_fx_NODU_hisp  <- sample(ind_NODU_hisp, size = round(length(ind_NODU_hisp) * NODU.fx["hisp"], 0))
      ppl_list[[t]]$fx[which(ppl_list[[t]]$ind %in% c(ind_fx_NODU_white, ind_fx_NODU_black, ind_fx_NODU_hisp))] <- 1
      m.tp <- trans.prob(ppl_list[[t]], params) # calculate the transition probabilities at cycle t
      n.nlx.mn$OEND <- n.nlx.mn$OEND * (1 - r.LossExp) + nx_avail_yr$OEND / 12
      n.nlx.mn$Pharm <- n.nlx.mn$Pharm * (1 - r.LossExp) + nx_avail_yr$Pharm / 12
    }

    set.seed((seed + t)*4) # set the seed for every individual for the random number generator
    ppl_list[[t]]$curr.state <- as.vector(samplev(probs = m.tp, m = 1)) # sample the next health state
    ind.oustate.chg <- filter(ppl_list[[t]], curr.state %in% v.oustate & OU.state != curr.state)$ind
    ppl_list[[t]]$OU.state[ind.oustate.chg] <- ppl_list[[t]]$curr.state[ind.oustate.chg]
    od_ppl <- ppl_list[[t]][ppl_list[[t]]$curr.state == "od", ]
    v.od[t] <- nrow(od_ppl)
    # ou.pop.resid <- ppl_list[[t]] %>% count(residence)
    crosstab_race_resid <- crosstab(ppl_list[[t]], row.vars = "residence", col.vars = "race", type = "f")$table
    crosstab_race_resid <- crosstab_race_resid[-nrow(crosstab_race_resid), -ncol(crosstab_race_resid)]
    crosstab_race_resid[crosstab_race_resid <= 0] <- 1 ## Replace population size with value 0 by 1
    crosstab_race_resid <- crosstab_race_resid[,c("white", "black", "hisp")]
    # Naloxone availability algorithm to determine the provability of naloxone used in a witnessed overdose
    p.nlx.avail.mx <- nlx.avail.algm(n.nlx = n.nlx.mn, crosstab_race_resid, OD_loc_pub, nlx.adj, cap, eff.pharNlx, strategy, t)
    decntree.out <- decision_tree(od_ppl, p.nlx.avail.mx, params, strategy, t, seed = seed + t)

    v.oddeath[t] <- sum(decntree.out[, "od.death"])
    v.oddeath.w[t] <- sum(decntree.out[decntree.out[ , "wtns"] == 1 , "od.death"])
    v.oddeath.w.race[t, "white"] <- sum(decntree.out[decntree.out[ , "wtns"] == 1 & decntree.out[,"race.ind"] == 1, "od.death"])
    v.oddeath.w.race[t, "black"] <- sum(decntree.out[decntree.out[ , "wtns"] == 1 & decntree.out[,"race.ind"] == 2, "od.death"])
    v.oddeath.w.race[t, "hisp"]  <- sum(decntree.out[decntree.out[ , "wtns"] == 1 & decntree.out[,"race.ind"] == 3, "od.death"])
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
    m.oddeath.fx.race[t,"white"] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$fx == 1 & od_ppl$race == "white", ])
    m.oddeath.fx.race[t,"black"] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$fx == 1 & od_ppl$race == "black", ])
    m.oddeath.fx.race[t,"hisp"]  <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$fx == 1 & od_ppl$race == "hisp", ])
    m.oddeath.op[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU", ])
    # m.oddeath.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state != "NODU" & od_ppl$OU.state != "preb", ])
    m.oddeath.st[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state == "NODU", ])
    m.EDvisits[t] <- sum(decntree.out[, "hospcare"])
    m.EDvisits.race[t, "white"] <- sum(decntree.out[decntree.out[,"race.ind"] == 1, "hospcare"])
    m.EDvisits.race[t, "black"] <- sum(decntree.out[decntree.out[,"race.ind"] == 2, "hospcare"])
    m.EDvisits.race[t, "hisp"]  <- sum(decntree.out[decntree.out[,"race.ind"] == 3, "hospcare"])
    
    ##ADDED for od deaths stratified by population groups##
    m.oddeath.preb[t]  <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state =="preb", ])
    m.oddeath.il.lr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state =="il.lr", ])
    m.oddeath.il.hr[t] <- nrow(od_ppl[od_ppl$curr.state == "dead" & od_ppl$OU.state =="il.hr", ])
    
    m.nlx.mn.OEND[t] <- sum(n.nlx.mn$OEND)
    m.nlx.mn.Pharm[t] <- sum(n.nlx.mn$Pharm)
    ####

    od.death.sum   <- od_ppl[od_ppl$curr.state == "dead", ] %>% count(residence)
    od.death.white <- od_ppl[od_ppl$curr.state == "dead" & od_ppl$race == "white", ] %>% count(residence)
    od.death.black <- od_ppl[od_ppl$curr.state == "dead" & od_ppl$race == "black", ] %>% count(residence)
    od.death.hisp  <- od_ppl[od_ppl$curr.state == "dead" & od_ppl$race == "hisp", ] %>% count(residence)
    for (dd in 1:nrow(od.death.sum)) {
      m.oddeath[t, od.death.sum$residence[dd]] <- od.death.sum$n[dd]
    }
    if(nrow(od.death.white) != 0){
      for (dd in 1:nrow(od.death.white)) {
        m.oddeath.white[t, od.death.white$residence[dd]] <- od.death.white$n[dd]
      }
    }
    if(nrow(od.death.black) != 0){
      for (dd in 1:nrow(od.death.black)) {
        m.oddeath.black[t, od.death.black$residence[dd]] <- od.death.black$n[dd]
      }
    }
    if(nrow(od.death.hisp) != 0){
      for (dd in 1:nrow(od.death.hisp)) {
        m.oddeath.hisp [t, od.death.hisp$residence[dd]]  <- od.death.hisp$n[dd]
      }
    }
    # ppl_list[[t]][which(ppl_list[[t]]$ind %in% od_ppl$ind), ] <- od_ppl
    # ppl_list[[t]][ppl_list[[t]]$curr.state == "od", ] <- od_ppl
    ppl_list[[t]][od_ppl$ind, ] <- od_ppl
    
    t0.eva <- (2020-2016 + 1)*12
    if (t > t0.eva){
      cost.matrix[t-t0.eva, ] <- Costs(state = ppl_list[[t]]$curr.state, OU.state = ppl_list[[t]]$OU.state, nlx = sum(unlist(nx_avail_yr)) / 12, count = list(n.EMS = n.EMS, n.hospcare = n.hospcare), params)
      qaly <- QALYs(state = ppl_list[[t]]$curr.state, OU.state = ppl_list[[t]]$OU.state, race = ppl_list[[t]]$race, params = params)
      qaly.matrix[t-t0.eva, "all"]   <- qaly$Q.TQ
      qaly.matrix[t-t0.eva, "white"] <- qaly$Q.white
      qaly.matrix[t-t0.eva, "black"] <- qaly$Q.black
      qaly.matrix[t-t0.eva, "hisp"] <- qaly$Q.hisp
    }
    ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] <- ppl_list[[t]]$age[ppl_list[[t]]$curr.state != "dead"] + ifelse(t %% 12 == 0, 1, 0) # update age for individuals that are still alive

    # record the individuals who die
    tmi <- rep(t, nrow(od_ppl[od_ppl$curr.state == "dead", ]))
    # dead_pop <- rbind(cbind(od_ppl[od_ppl$curr.state == "dead", ], tmi), dead_pop)

    # replace deceased individuals with ones with the same initial characteristics (ever.od reset as 0), disabled for MA model, only during calibration period
    if (t <= t0.eva){
      ppl_list[[t]]$age[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.age[ppl_list[[t]]$curr.state == "dead"]
      ppl_list[[t]]$ever.od[ppl_list[[t]]$curr.state == "dead"] <- 0
      ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state != "inact"]
      ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"] <- ppl_list[[t]]$OU.state[ppl_list[[t]]$curr.state == "dead" & ppl_list[[t]]$init.state == "inact"]
      ppl_list[[t]]$curr.state[ppl_list[[t]]$curr.state == "dead"] <- ppl_list[[t]]$init.state[ppl_list[[t]]$curr.state == "dead"]
    }
    # # Allow fixed number of new population entering the model (to overcome high-risk groups dying too fast) without affecting QALY and cost estimates
    # ppl_list[[t]] <- rbind(ppl_list[[t]], add_pop[((t-1)*n.add.pop+1):(t*n.add.pop), ])
    
    # cat('\r', paste(round(t/timesteps * 100), "% done", sep = " "))       # display the progress of the simulation
  } # end the loop for the time steps


  total.cost <- sum(cost.matrix[, "TotalCost"] * v.dwc) # total (discounted) cost
  total.qaly <- sum(qaly.matrix[, "all"] * v.dwc) # total (discounted) qalys among all population
  total.qaly.w <- sum(qaly.matrix[, "white"] * v.dwc) # total (discounted) qalys among white
  total.qaly.b <- sum(qaly.matrix[, "black"] * v.dwc) # total (discounted) qalys among black
  total.qaly.h <- sum(qaly.matrix[, "hisp"] * v.dwc) # total (discounted) qalys among hispanic

  if (PT.out == TRUE) {
    pop.trace <- ppl_list
  } else {
    pop.trace <- NULL
  }
  # if (grepl("10K", strategy)){
  #   n.nlx.OEND.str.base <- n.nlx.str$OEND[, -ncol(n.nlx.str$OEND)]
  #   n.nlx.OEND.str.add.total <- n.nlx.str$OEND[, ncol(n.nlx.str$OEND)]
  #   if (strategy == "SSP_10K"){
  #     SSP_drug_resid <- crosstab_drug_resid[, c("il.hr", "NODU")]*rep(c(1, 0.133), each = nrow(crosstab_drug_resid))
  #     n.nlx.OEND.str <- round(rowSums(n.nlx.OEND.str.add.total * SSP_drug_resid / sum(SSP_drug_resid)), 0) + n.nlx.OEND.str.base
  #   # } else if (strategy == "StreetOutreach_10K"){
  #   #   StreetOutreach_drug_resid <- crosstab_drug_resid[, c("il.lr", "il.hr", "NODU")]
  #   #   n.nlx.OEND.str <- round(rowSums(n.nlx.OEND.str.add.total * StreetOutreach_drug_resid / sum(StreetOutreach_drug_resid)), 0) + n.nlx.OEND.str.base
  #   } else if (strategy == "MailEvent_10K"){
  #     MailEvent_drug_resid <- crosstab_drug_resid
  #     n.nlx.OEND.str <- round(rowSums(n.nlx.OEND.str.add.total * MailEvent_drug_resid / sum(MailEvent_drug_resid)), 0) + n.nlx.OEND.str.base
  #   } else if (strategy == "Healthcare_10K"){
  #     Healthcare_drug_resid <- crosstab_drug_resid[, c("preb")]
  #     n.nlx.OEND.str <- round(n.nlx.OEND.str.add.total * Healthcare_drug_resid / sum(Healthcare_drug_resid), 0) + n.nlx.OEND.str.base
  #   }
  # } else {
  #   n.nlx.OEND.str <- n.nlx.str$OEND
  # }
  n.nlx.OEND.str <- n.nlx.str$OEND
  print("Saving results")
  results <- list(
    # dead_pop = dead_pop,
    v.oddeath = v.oddeath, m.oddeath = m.oddeath, m.oddeath.white = m.oddeath.white, m.oddeath.black = m.oddeath.black, m.oddeath.hisp = m.oddeath.hisp, v.od = v.od,
    cost.matrix = cost.matrix, total.cost = total.cost,
    qaly.matrix = qaly.matrix, total.qaly = total.qaly, total.qaly.w = total.qaly.w, total.qaly.b = total.qaly.b, total.qaly.h = total.qaly.h, 
    pop.trace = pop.trace,  n.nlx.OEND.str = n.nlx.OEND.str,
    m.oddeath.fx = m.oddeath.fx, m.oddeath.fx.race = m.oddeath.fx.race,
    m.oddeath.op = m.oddeath.op, m.oddeath.st = m.oddeath.st, m.oddeath.hr = m.oddeath.hr, 
    m.EDvisits = m.EDvisits, m.EDvisits.race = m.EDvisits.race,
    v.odpriv = v.odpriv, v.odpubl = v.odpubl, v.deathpriv = v.deathpriv, v.deathpubl = v.deathpubl, v.nlxused = v.nlxused, v.oddeath.w = v.oddeath.w, v.oddeath.w.race = v.oddeath.w.race,
    m.oddeath.preb = m.oddeath.preb, m.oddeath.il.lr = m.oddeath.il.lr, m.oddeath.il.hr = m.oddeath.il.hr, m.nlx.mn.OEND = m.nlx.mn.OEND, m.nlx.mn.Pharm = m.nlx.mn.Pharm, 
    end.pop.state = table(ppl_list[[t]]$curr.state)
  ) # store the results from the simulation in a list
  return(results) # return the results
} # end of the MicroSim function
