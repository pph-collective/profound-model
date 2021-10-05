#' Create dataframes for params and model outputs
#'
#' @description
#' `input_setup()` adds necessary parameters to the model's input params
#'
#' `output_setup()` creates a dataframe to store the model's output data
#'
#' @param params The parameters for the model.
#' @param data Empirical data to inform the model.
#'
#' @returns
#' `input_setup()` returns a named list of model parameters
#'
#' `output_setup()` returns a dataframe to store model outputs
#'

library("yaml")
input_setup <- function(params, data) {
  params$pop.info <- c(
    "sex", "race", "age", "residence", "curr.state",
    "OU.state", "init.age", "init.state", "ever_od", "fx"
  ) # information for each model individual
  params$agent_states <- c(
    "rx", "il_lr", "il_hr", "inact", "NODU", "relap", "dead"
  )
  params$oustate <- c("rx", "il_lr", "il_hr") # opioid use states
  params$num_years <- params$year_end - params$year_start + 1
  params$timesteps <- 12 * params$num_years # number of time cycles (in month)
  params$regions <- data$regions
  return(params)
}

output_setup <- function(params) {
  output <- data.frame(
    t = 1:params$timesteps,
    v.od = rep(0, times = params$timesteps) # count of overdoses at timestep
  )

  output$oddeath <- NA # count of overdose deaths
  output$death_wtns <- NA # count of witnessed deaths
  output$total_cost <- NA
  output$nx_cost <- NA
  output$priv_od <- NA # count of overdose events at private setting
  output$vpub_od <- NA # count of overdose events at public setting
  output$priv_death <- NA # count of overdose deaths at private setting
  output$pub_death <- NA # count of overdose deaths at public setting
  output$nlxused <- NA # count of naloxone kits used
  output$fx_deaths <- NA # count of overdose deaths with fentanyl present
  output$oud_deaths <- NA # count of overdose deaths among opioid users
  output$stim_deaths <- NA # count of overdose deaths among stimulant users
  output$edvisits <- NA # count of opioid overdose-related ED visits
  output$hr_deaths <- NA # count of overdose deaths among high-risk opioid

  for (region in params$regions) {
    output[region] <- NA
  }

  return(output)
}
