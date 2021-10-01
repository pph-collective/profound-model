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
    "ever_od", "fx")

  # get actual population size
  ppl_size <- round(ppl_size * (oud_prev + mean(c(nodu_m_prev, nodu_f_prev))))
  init_ppl <- data.frame(matrix(0, (ppl_size), length(ppl_info)))
  colnames(init_ppl) <- ppl_info

  for (i in 1:ppl_size) {
    set.seed(seed + i)
    agent <- list()
    agent$residence <- sample(p_region$region, 1, prob = p_region$prob)
    # get age/race/sex
    age_probs <- p_region[agent$residence, -c(1, 2)]
    dem <- strsplit(sample(colnames(age_probs), 1, prob = age_probs), "_")
    agent$sex <- dem[[1]][1]
    agent$race <- dem[[1]][2]
    agent$age <- dem[[1]][3]

    oud_prob <- oud_demo$pe[
      oud_demo$age == agent$age &
      oud_demo$sex == agent$sex & oud_demo$race == agent$race]

    stim_prob <- stim_demo$pe[
      stim_demo$sex == agent$sex &
      stim_demo$race == agent$race & stim_demo$age == agent$age]

    oud_prob <- oud_prob / (oud_prob + stim_prob)

    agent$ever_od <- FALSE
    # OUD agent?
    if (runif(1) < oud_prob) {
      agent <- make_oud_agent(agent, initials, agent_states)
    } else {
      agent$curr.state <- agent$init.state <- agent$OU.state <- "NODU"
    }

    age <- agent$age
    agent$age <- round(runif(1, as.integer(substr(age, 1, 2)), as.integer(substr(age, nchar(age) - 1, nchar(age)))))
    agent$init.age <- agent$age
    agent$fx <- 0  # TO_REVIEW seems like this is how it is in the original?
    # add agent to pop
    for (j in names(agent)) {
      init_ppl[i, j] <- agent[j]
    }
  }
  return(init_ppl)
}

make_oud_agent <- function(agent, params, agent_states) {
  list2env(params, environment())
  sex <- agent$sex
  oud.state <- agent_states[1:4]
  # non-illicit drug use probability
  non_il <- (1 - init_inactive) * (1 - get(paste0("init_il_", sex)))

  # high risk use probability
  hr <- (1 - init_inactive) * get(paste0("init_il_", sex)) * get(paste0("init_il_hr_", sex))

  # illicit, non-high risk probability
  il <- (1 - init_inactive) * get(paste0("init_il_", sex)) * (1 - get(paste0("init_il_hr_", sex)))

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
  p_ever_od <- get(paste0("init_everod_", agent$OU.state))
  agent$ever_od <- runif(1) < p_ever_od
  return(agent)
}
