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

initiate_ppl <- function(data, agent_states, seed = 2021) {
  initials <- data$initials
  list2env(initials, environment())
  ppl_info <- c(
    "sex", "race", "age",
    "residence",
    "curr.state", "OU.state",
    "init.age", "init.state",
    "ever.od", "fx")

  # get actual population size
  ppl_size <- round(ppl_size * (prev.oud + mean(c(prev.NODU.m, prev.NODU.f))))
  init_ppl <- data.frame(matrix(0, (ppl_size), length(ppl_info)))
  colnames(init_ppl) <- ppl_info

  print(ppl_size)

  for (i in 1:ppl_size) {
    tic(paste0(i, ": "))
    set.seed(seed + i)
    agent <- list()
    agent$residence <- sample(p_region$region, 1, prob = p_region$prob)

    # get age/race/sex
    age_probs <- p_region[agent$residence, -c(1, 2)]
    dem <- strsplit(sample(colnames(age_probs), 1, prob = age_probs), "_")
    agent$sex <- dem[[1]][1]
    agent$race <- dem[[1]][2]
    agent$age <- dem[[1]][3]
    agent$init.age <- agent$age

    oud_prob <- OUDDemo$pe[
      OUDDemo$age == agent$age &
      OUDDemo$sex == agent$sex & OUDDemo$race == agent$race]

    stim_prob <- StimDemo$pe[
      StimDemo$sex == agent$sex &
      StimDemo$race == agent$race & StimDemo$age == agent$age]

    oud_prob <- oud_prob / (oud_prob + stim_prob)

    agent$ever.od <- FALSE
    # OUD agent?
    if (runif(1) < oud_prob) {
      agent <- make_oud_agent(agent, initials, agent_states)
    } else {
      agent$curr.state <- agent$init.state <- agent$OU.state <- "NODU"
    }
    age <- agent$age
    agent$age <- round(runif(1, as.integer(substr(age, 1, 2)), as.integer(substr(age, nchar(age) - 1, nchar(age)))))

    agent$fx <- 0  # TO_REVIEW seems like this is how it is in the original?
    # add agent to pop
    rbind(init_ppl, agent)
    toc()
  }

  return(init_ppl)
}

make_oud_agent <- function(agent, params, agent_states) {
  list2env(params, environment())
  sex <- agent$sex
  oud.state <- agent_states[1:4]
  # non-illicit drug use probability
  non_il <- (1 - init_inactive) * (1 - get(paste0("ini.il.", sex)))

  # high risk use probability
  hr <- (1 - init_inactive) * get(paste0("ini.il.", sex)) * get(paste0("ini.il.hr.", sex))

  # illicit, non-high risk probability
  il <- (1 - init_inactive) * get(paste0("ini.il.", sex)) * (1 - get(paste0("ini.il.hr.", sex)))

  # create probability vector
  oud.prob <- c(non_il, il, hr, init_inactive)

  # sample to find current state
  agent$curr.state <- sample(oud.state, size = 1, prob = oud.prob)

  agent$init.state <- agent$curr.state

  agent$OU.state <- agent$curr.state


  # check what type of drugs used if currently inactive
  if (agent$curr.state == "inact") {
    agent$OU.state <- sample(oud.state[1:3], size = 1, prob = oud.prob[1:3])
  }

  # has the agent previously overdosed?
  p_ever_od <- get(paste0("ini.everod.", agent$OU.state))
  agent$ever.od <- runif(1) < p_ever_od

  return(agent)
}
