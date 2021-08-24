###########################################################################################
#########################       Input and Output setup     ################################
###########################################################################################

###########################################################################################
# This is to define all the input and output variables required in the model (consistent across all scenarios)

library("yaml")
# INPUT setup
input_setup <- function(params, data){
  params$pop.info <- c(
    "sex", "race", "age", "residence", "curr.state",
    "OU.state", "init.age", "init.state", "ever.od", "fx"
  ) # information for each model individual
  params$agent_states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead") # vector for state names
  params$v.oustate <- c("preb", "il.lr", "il.hr") # vector for active opioid use state names
  num_states <- length(params$agent_states) # number of states
  params$num_years <- params$year_end - params$year_start + 1
  params$timesteps <- 12 * params$num_years # number of time cycles (in month)
  params$v.region <- data$v.region
  num_regions <- length(params$v.region) # number of regions
  params$strategies = c("SQ", "expand", "program")
  return(params)
}

# OUTPUT matrices and vectors
output_setup <- function(params){
  cost_labels = c("TotalCost", "NxCost")
  output <- data.frame(
    v.od = rep(0, times = params$timesteps) # count of overdoses at timestep
  )

  output$v.oddeath = NA # count of overdose deaths at each time step
  output$v.oddeath.w = NA # count of overdose deaths that were witnessed at each time step
  output$total_cost = NA
  output$nx_cost = NA
  output$v.odpriv = NA # count of overdose events occurred at private setting at each time step
  output$v.odpubl = NA # count of overdose events occurred at public setting at each time step
  output$v.deathpriv = NA # count of overdose deaths occurred at private setting at each time step
  output$v.deathpubl = NA # count of overdose deaths occurred at public setting at each time step
  output$v.nlxused = NA # count of naloxone kits used at each time step
  output$m.oddeath.fx = NA # count of overdose deaths with fentanyl present at each time step
  output$m.oddeath.op = NA # count of overdose deaths among opioid users at each time step
  output$m.oddeath.st = NA # count of overdose deaths among stimulant users at each time step
  output$m.EDvisits = NA # count of opioid overdose-related ED visits at each time step
  output$m.oddeath.hr = NA # count of overdose deaths among high-risk opioid users (inject heroin) at each time step

  for (region in params$v.region) {
    output[region] <- NA
  }


  return(output)
}

