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
#################        Transistion probability function      #########################
########################################################################################
library(dplyr)
# TODO there must be a way to do all of the filters at once?
# the issue seems to be the age -- might want to change this to allow age category as well

trans.prob <- function(pop.t, params, data) {
  # initialize a matrix to store transition probabilities of each individual
  trans.prob.matrix <- matrix(NA, length(params$agent_states) + 1, nrow(pop.t)) # add one row for OD
  rownames(trans.prob.matrix) <- c(params$agent_states, "od")

  # create a vector to store baseline mortality (excluding od death) for each individual according to age and treatment
  mor.rate <- numeric(nrow(pop.t))

  tic("mor rate 1")
  mor.rate[filter(pop.t, age %in% c(10:14))$ind] <- data$mortality_probs["drug", "10to14"]
  mor.rate[filter(pop.t, age %in% c(15:24))$ind] <- data$mortality_probs["drug", "15to24"]
  mor.rate[filter(pop.t, age %in% c(25:34))$ind] <- data$mortality_probs["drug", "25to34"]
  mor.rate[filter(pop.t, age %in% c(35:44))$ind] <- data$mortality_probs["drug", "35to44"]
  mor.rate[filter(pop.t, age %in% c(45:54))$ind] <- data$mortality_probs["drug", "45to54"]
  mor.rate[filter(pop.t, age %in% c(55:64))$ind] <- data$mortality_probs["drug", "55to64"]
  mor.rate[filter(pop.t, age >= 65)$ind] <- data$mortality_probs["drug", "65over"]
  toc()

  # create a vector to store probability of overdose for each individual according to ever overdosed and fentanyl
  # TO_REVIEW what does "multi" mean here
  od.rate <- numeric(nrow(pop.t))
  od.rate[filter(pop.t, curr.state == "rx" & ever_od == 0 & fx == 0)$ind] <- data$overdose_probs["rx", "first"]
  od.rate[filter(pop.t, curr.state == "rx" & ever_od == 1 & fx == 0)$ind] <- data$overdose_probs["rx", "subs"]
  od.rate[filter(pop.t, curr.state == "rx" & ever_od == 0 & fx == 1)$ind] <- data$overdose_probs["rx", "first"] * data$multi_fx
  od.rate[filter(pop.t, curr.state == "rx" & ever_od == 1 & fx == 1)$ind] <- data$overdose_probs["rx", "subs"] * data$multi_fx
  od.rate[filter(pop.t, curr.state == "il_lr" & ever_od == 0 & fx == 0)$ind] <- data$overdose_probs["il_lr", "first"]
  od.rate[filter(pop.t, curr.state == "il_lr" & ever_od == 1 & fx == 0)$ind] <- data$overdose_probs["il_lr", "subs"]
  od.rate[filter(pop.t, curr.state == "il_lr" & ever_od == 0 & fx == 1)$ind] <- data$overdose_probs["il_lr", "first"] * data$multi_fx
  od.rate[filter(pop.t, curr.state == "il_lr" & ever_od == 1 & fx == 1)$ind] <- data$overdose_probs["il_lr", "subs"] * data$multi_fx
  od.rate[filter(pop.t, curr.state == "il_hr" & ever_od == 0 & fx == 0)$ind] <- data$overdose_probs["il_hr", "first"]
  od.rate[filter(pop.t, curr.state == "il_hr" & ever_od == 1 & fx == 0)$ind] <- data$overdose_probs["il_hr", "subs"]
  od.rate[filter(pop.t, curr.state == "il_hr" & ever_od == 0 & fx == 1)$ind] <- data$overdose_probs["il_hr", "first"] * data$multi_fx
  od.rate[filter(pop.t, curr.state == "il_hr" & ever_od == 1 & fx == 1)$ind] <- data$overdose_probs["il_hr", "subs"] * data$multi_fx
  od.rate[filter(pop.t, curr.state == "NODU" & ever_od == 0 & fx == 0)$ind] <- data$overdose_probs["NODU", "first"]
  od.rate[filter(pop.t, curr.state == "NODU" & ever_od == 1 & fx == 0)$ind] <- data$overdose_probs["NODU", "subs"]
  od.rate[filter(pop.t, curr.state == "NODU" & ever_od == 0 & fx == 1)$ind] <- data$overdose_probs["il_lr", "first"] * data$multi_NODU_fx * data$multi_fx
  od.rate[filter(pop.t, curr.state == "NODU" & ever_od == 1 & fx == 1)$ind] <- data$overdose_probs["il_lr", "subs"] * data$multi_NODU_fx * data$multi_fx
  od.rate[filter(pop.t, curr.state == "relap" & ever_od == 0 & OU.state == "rx")$ind] <- data$overdose_probs["rx", "first"] * data$multi_relap
  od.rate[filter(pop.t, curr.state == "relap" & ever_od == 1 & OU.state == "rx")$ind] <- data$overdose_probs["rx", "subs"] * data$multi_relap
  od.rate[filter(pop.t, curr.state == "relap" & ever_od == 0 & OU.state == "il_lr")$ind] <- data$overdose_probs["il_lr", "first"] * data$multi_relap
  od.rate[filter(pop.t, curr.state == "relap" & ever_od == 1 & OU.state == "il_hr")$ind] <- data$overdose_probs["il_lr", "subs"] * data$multi_relap
  od.rate[filter(pop.t, curr.state == "relap" & ever_od == 0 & OU.state == "il_hr")$ind] <- data$overdose_probs["il_hr", "first"] * data$multi_relap
  od.rate[filter(pop.t, curr.state == "relap" & ever_od == 1 & OU.state == "il_hr")$ind] <- data$overdose_probs["il_hr", "subs"] * data$multi_relap
  # update the trans.prob matrix with the corresponding probabilities
  ind.rx <- pop.t$curr.state == "rx"

  if (sum(ind.rx) != 0) {
    trans.prob.matrix[, ind.rx] <- rbind(
      1 - data$p_rx2il_lr - data$p_rx2inact - mor.rate[ind.rx] - od.rate[ind.rx],
      rep(data$p_rx2il_lr, sum(ind.rx)),
      rep(0, sum(ind.rx)),
      rep(data$p_rx2inact, sum(ind.rx)),
      rep(0, sum(ind.rx)),
      rep(0, sum(ind.rx)),
      mor.rate[ind.rx],
      od.rate[ind.rx]
    )
  }

  ind.il_lr <- pop.t$curr.state == "il_lr"
  if (sum(ind.il_lr) != 0) {
    trans.prob.matrix[, ind.il_lr] <- rbind(
      rep(0, sum(ind.il_lr)),
      1 - data$p_il_lr2il_hr - data$p_il_lr2inact - mor.rate[ind.il_lr] - od.rate[ind.il_lr],
      rep(data$p_il_lr2il_hr, sum(ind.il_lr)),
      rep(data$p_il_lr2inact, sum(ind.il_lr)),
      rep(0, sum(ind.il_lr)),
      rep(0, sum(ind.il_lr)),
      mor.rate[ind.il_lr],
      od.rate[ind.il_lr]
    )
  }

  ind.il_hr <- pop.t$curr.state == "il_hr"
  if (sum(ind.il_hr) != 0) {
    trans.prob.matrix[, ind.il_hr] <- rbind(
      rep(0, sum(ind.il_hr)),
      rep(data$p_il_hr2il_lr, sum(ind.il_hr)),
      1 - data$p_il_hr2il_lr - data$p_il_hr2il_lr - mor.rate[ind.il_hr] - od.rate[ind.il_hr],
      rep(data$p_il_hr2il_lr, sum(ind.il_hr)),
      rep(0, sum(ind.il_hr)),
      rep(0, sum(ind.il_hr)),
      mor.rate[ind.il_hr],
      od.rate[ind.il_hr]
    )
  }

  ind.inact <- pop.t$curr.state == "inact"
  if (sum(ind.inact) != 0) {
    trans.prob.matrix[, ind.inact] <- rbind(
      rep(0, sum(ind.inact)),
      rep(0, sum(ind.inact)),
      rep(0, sum(ind.inact)),
      1 - data$p_inact2relap - mor.rate[ind.inact] - od.rate[ind.inact],
      rep(0, sum(ind.inact)),
      rep(data$p_inact2relap, sum(ind.inact)),
      mor.rate[ind.inact],
      od.rate[ind.inact]
    )
  }

  ind.NODU <- pop.t$curr.state == "NODU"
  if (sum(ind.NODU) != 0) {
    trans.prob.matrix[, ind.NODU] <- rbind(
      rep(0, sum(ind.NODU)),
      rep(0, sum(ind.NODU)),
      rep(0, sum(ind.NODU)),
      rep(0, sum(ind.NODU)),
      1 - mor.rate[ind.NODU] - od.rate[ind.NODU],
      rep(0, sum(ind.NODU)),
      mor.rate[ind.NODU],
      od.rate[ind.NODU]
    )
  }

  ind.relap <- pop.t$curr.state == "relap"
  if (sum(ind.relap) != 0) {
    OU.v <- filter(pop.t, curr.state == "relap")$OU.state
    num_states <- length(params$agent_states)
    relap.m <- matrix(0, num_states + 1, sum(ind.relap))
    for (r in 1:sum(ind.relap)) {
      relap.m[which(OU.v[r] == params$agent_states), r] <- 1 - mor.rate[ind.relap][r] - od.rate[ind.relap][r]
      relap.m[which("dead" == params$agent_states), r] <- mor.rate[ind.relap][r]
      relap.m[num_states + 1, r] <- od.rate[ind.relap][r]
    }
    trans.prob.matrix[, ind.relap] <- relap.m
  }

  trans.prob.matrix[, pop.t$curr.state == "dead"] <- c(0, 0, 0, 0, 0, 0, 1, 0)

  return(t(trans.prob.matrix))
}

# random sampling for vector of outcomes
# TODO rewrite to store stats
samplev <- function(probs, m) {
  # get dimensions of probabilities
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  states <- dimnames(probs)[[2]]

  if (!length(states)) {
    states <- 1:k
  }
  ran <- matrix(states[1], ncol = m, nrow = n)

  U <- t(probs)

  for (i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }

  if (any((U[k, ] - 1) > 1e-05)) {
    stop("error in multinom: probabilities do not sum to 1")
  }

  for (j in 1:m) {
    # TO_REVIEW what is un?
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- states[1 + colSums(un > U)]
  }

  return(ran)
}
