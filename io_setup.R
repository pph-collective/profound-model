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
  output <- data.frame(
    v.od = rep(0, times = params$timesteps) # count of overdoses at timestep
  )

  output$v.oddeath = NA # count of overdose deaths
  output$v.oddeath.w = NA # count of witnessed deaths
  output$total_cost = NA
  output$nx_cost = NA
  output$v.odpriv = NA # count of overdose events at private setting
  output$v.odpubl = NA # count of overdose events at public setting
  output$v.deathpriv = NA # count of overdose deaths at private setting
  output$v.deathpubl = NA # count of overdose deaths at public setting
  output$v.nlxused = NA # count of naloxone kits used
  output$m.oddeath.fx = NA # count of overdose deaths with fentanyl present
  output$m.oddeath.op = NA # count of overdose deaths among opioid users
  output$m.oddeath.st = NA # count of overdose deaths among stimulant users
  output$m.EDvisits = NA # count of opioid overdose-related ED visits
  output$m.oddeath.hr = NA # count of overdose deaths among high-risk opioid users (inject heroin)

  for (region in params$v.region) {
    output[region] <- NA
  }

  return(output)
}

