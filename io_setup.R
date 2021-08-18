###########################################################################################
#########################       Input and Output setup     ################################
###########################################################################################

###########################################################################################
# This is to define all the input and output variables required in the model (consistent across all scenarios)

library("yaml")
# INPUT setup
input_setup <- function(params){
  params$pop.info <- c(
    "sex", "race", "age", "residence", "curr.state",
    "OU.state", "init.age", "init.state", "ever.od", "fx"
  ) # information for each model individual
  params$agent_states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead") # vector for state names
  params$v.oustate <- c("preb", "il.lr", "il.hr") # vector for active opioid use state names
  num_states <- length(params$agent_states) # number of states
  params$num_years <- params$year_end - params$year_start + 1
  params$timesteps <- 12 * (params$year_end - params$year_start) # number of time cycles (in month)
  num_regions <- length(params$v.region) # number of regions
  return(params)
}

# OUTPUT matrices and vectors
output_setup <- function(params){
  output <- list()
  output$v.od <- rep(0, times = params$timesteps) # count of overdose events at each time step
  output$v.oddeath <- rep(0, times = params$timesteps) # count of overdose deaths at each time step
  output$v.oddeath.w <- rep(0, times = params$timesteps) # count of overdose deaths that were witnessed at each time step
  output$m.oddeath <- matrix(0, nrow = params$timesteps, ncol = length(params$v.region))
  colnames(output$m.oddeath) <- params$v.region
  output$v.odpriv <- rep(0, times = params$timesteps) # count of overdose events occurred at private setting at each time step
  output$v.odpubl <- rep(0, times = params$timesteps) # count of overdose events occurred at public setting at each time step
  output$v.deathpriv <- rep(0, times = params$timesteps) # count of overdose deaths occurred at private setting at each time step
  output$v.deathpubl <- rep(0, times = params$timesteps) # count of overdose deaths occurred at public setting at each time step
  output$v.nlxused <- rep(0, times = params$timesteps) # count of naloxone kits used at each time step
  output$strategies <- c("SQ", "expand", "program") # store the strategy names
  cost_labels <- c("TotalCost", "NxCost")
  output$cost.matrix <- matrix(0, nrow = params$timesteps, ncol = length(cost_labels))
  colnames(output$cost.matrix) <- cost_labels
  output$m.oddeath.fx <- rep(0, times = params$timesteps) # count of overdose deaths with fentanyl present at each time step
  output$m.oddeath.op <- rep(0, times = params$timesteps) # count of overdose deaths among opioid users at each time step
  output$m.oddeath.st <- rep(0, times = params$timesteps) # count of overdose deaths among stimulant users at each time step
  output$m.EDvisits <- rep(0, times = params$timesteps) # count of opioid overdose-related ED visits at each time step
  output$m.oddeath.hr <- rep(0, times = params$timesteps) # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
  return(output)
}
