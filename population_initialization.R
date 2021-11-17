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
#################        Population creation/initialization function      ##############
# To initialize study population at time = 0, define all individual attributes
########################################################################################

initiate_ppl <- function(inits, agent_states, seed = 2021) {
  ppl_info <- c(
    "sex", "race", "age",
    "residence",
    "curr.state", "OU.state",
    "init.age", "init.state",
    "ever.od", "fx")
  ## opioid population
  # Determine people with opioid use disorder as estimated by NSDUH
  oud_pop_nsduh <- inits$demo.mx * inits$oud_demo # number of OUD population estimates for each subgroup in each region according to NSDUH prevalence estimates
  oud.adj <-  inits$oud_prev / (sum(oud_pop_nsduh) / sum(inits$demo.mx)) # adjuster for OUD population according to overall oud prevalence and estimated one from the previous step
  oud_pop <- oud_pop_nsduh * oud.adj
  prop.oud.region <- colSums(oud_pop) / sum(oud_pop)
  oud_demo <- oud_pop / matrix(rep(colSums(oud_pop), each = nrow(oud_pop)), nrow = nrow(oud_pop), ncol = ncol(oud_pop))

  ## non-opioid population
  stim_pop <- inits$demo.mx * inits$stim_demo # number of NODU population estimates for each subgroup in each region according to NSDUH prevalence estimates
  prop.nodu.region <- colSums(stim_pop) / sum(stim_pop)
  nodu_demo <- stim_pop / matrix(rep(colSums(stim_pop), each = nrow(stim_pop)), nrow = nrow(stim_pop), ncol = ncol(stim_pop))

  ## initialize the population matrix
  total.oud <- round(sum(oud_pop), 0)
  total.nodu <- round(sum(stim_pop), 0)
  init_ppl <- data.frame(matrix(0, (total.oud + total.nodu), length(ppl_info) + 1)) # initial population matrix
  dimnames(init_ppl)[[2]] <- c("ind", ppl_info)
  init_ppl$ind <- 1:(total.oud + total.nodu)
  set.seed(seed)
  oud.region.ind <- sample(1:length(inits$regions), size = total.oud, prob = prop.oud.region, replace = TRUE)
  nodu.region.ind <- sample(1:length(inits$regions), size = total.nodu, prob = prop.nodu.region, replace = TRUE)

  ## OPIOID USE POPULATION
  for (i in 1:total.oud) {
    set.seed(seed + i)
    # determine demographic
    demo.ind <- sample(1:nrow(oud_demo), size = 1, prob = oud_demo[, oud.region.ind[i]])
    sex <- inits$demographic$sex[demo.ind]
    race <- inits$demographic$race[demo.ind]
    age.ind <- inits$demographic$age[demo.ind]
    
    if (age.ind == "12to17") {
      age <- init.age <- sample(12:17, size = 1)
    } else if (age.ind == "18to25") {
      age <- init.age <- sample(18:25, size = 1)
    } else if (age.ind == "26to34") {
      age <- init.age <- sample(26:34, size = 1)
    } else if (age.ind == "35to49") {
      age <- init.age <- sample(35:49, size = 1)
    } else if (age.ind == "50to64") {
      age <- init.age <- sample(50:64, size = 1)
    } else if (age.ind == "65older") {
      age <- init.age <- 65
    }

    # determine opioid use state
    oud.state <- agent_states[1:4]
    inact <- inits$init_inactive
    # male probability
    il <- inits$init_il_m
    hr <- inits$init_il_hr_m
    oud_prob_m <- c((1 - inact) * (1 - il), (1 - inact) * il * (1 - hr), (1 - inact) * il * hr, inact)
    # female probability
    il <- inits$init_il_f
    hr <- inits$init_il_hr_f
    oud_prob_f <- c((1 - inact) * (1 - il), (1 - inact) * il * (1 - hr), (1 - inact) * il * hr, inact)
    if (sex == "m") {
      curr.state <- init.state <- oud.state[sample(1:length(oud.state), size = 1, prob = oud_prob_m)]
    } else {
      curr.state <- init.state <- oud.state[sample(1:length(oud.state), size = 1, prob = oud_prob_f)]
    }

    if (curr.state == "inact") {
      if (sex == "m") {
        OU.state <- oud.state[sample(1:length(oud.state[1:3]), size = 1, prob = oud_prob_m[1:3] / sum(oud_prob_m[1:3]))]
      } else {
        OU.state <- oud.state[sample(1:length(oud.state[1:3]), size = 1, prob = oud_prob_f[1:3] / sum(oud_prob_f[1:3]))]
      }
    } else {
      OU.state <- curr.state
    }

    # # determine fentanyl use
    # fx         <- sample(0:1, size = 1, prob = c(1-ini.OUD.fx, ini.OUD.fx))

    # determine ever overdosed
    if (OU.state == "rx") {
      ever.od <- sample(0:1, size = 1, prob = c(1 - inits$init_everod_rx, inits$init_everod_rx))
    } else if (OU.state == "il_lr") {
      ever.od <- sample(0:1, size = 1, prob = c(1 - inits$init_everod_il_lr, inits$init_everod_il_lr))
    } else if (OU.state == "il_hr") {
      ever.od <- sample(0:1, size = 1, prob = c(1 - inits$init_everod_il_hr, inits$init_everod_il_hr))
    } else {
      print(OU.state)
      stopifnot(FALSE)
    }


    init_ppl$sex[i] <- sex
    init_ppl$race[i] <- race
    init_ppl$age[i] <- age
    init_ppl$residence[i] <- inits$regions[oud.region.ind[i]]
    init_ppl$curr.state[i] <- curr.state
    init_ppl$OU.state[i] <- OU.state
    init_ppl$init.age[i] <- init.age
    init_ppl$init.state[i] <- init.state
    init_ppl$ever.od[i] <- ever.od
    # init_ppl$fx[i]         <- fx
  }

  ## NON-OPIOID USE POPULATION
  for (i in 1:total.nodu) {
    set.seed(seed + i)

    curr.state <- init.state <- OU.state <- "NODU"

    # determine demographic
    demo.ind <- sample(1:nrow(nodu_demo), size = 1, prob = nodu_demo[, nodu.region.ind[i]])
    # sex <- v.demo.sex[demo.ind]
    sex <- inits$demographic$sex[demo.ind]
    race <- inits$demographic$race[demo.ind]
    age.ind <- inits$demographic$age[demo.ind]
    if (age.ind == "12to17") {
      age <- init.age <- sample(12:17, size = 1)
    } else if (age.ind == "18to25") {
      age <- init.age <- sample(18:25, size = 1)
    } else if (age.ind == "26to34") {
      age <- init.age <- sample(26:34, size = 1)
    } else if (age.ind == "35to49") {
      age <- init.age <- sample(35:49, size = 1)
    } else if (age.ind == "50to64") {
      age <- init.age <- sample(50:64, size = 1)
    } else if (age.ind == "65older") {
      age <- init.age <- 65
    }

    # # determine fentanyl use
    # fx         <- sample(0:1, size = 1, prob = c(1-ini.NOUD.fx, ini.NOUD.fx))

    # determine ever overdosed
    ever.od <- sample(0:1, size = 1, prob = c(1 - inits$init_everod_sti, inits$init_everod_sti))


    init_ppl$sex[i + total.oud] <- sex
    init_ppl$race[i + total.oud] <- race
    init_ppl$age[i + total.oud] <- age
    init_ppl$residence[i + total.oud] <- inits$regions[nodu.region.ind[i]]
    init_ppl$curr.state[i + total.oud] <- curr.state
    init_ppl$OU.state[i + total.oud] <- OU.state
    init_ppl$init.age[i + total.oud] <- init.age
    init_ppl$init.state[i + total.oud] <- init.state
    init_ppl$ever.od[i + total.oud] <- ever.od
    # init_ppl$fx[i+total.oud]         <- fx
  }
  return(init_ppl)
}
