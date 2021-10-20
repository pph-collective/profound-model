# TODO file header

parallel.fun <- function(calib.seed, params) {
  ## Run the simulation model over all sets of parameters
  source("io_setup.R")
  output <- output_setup(params)
  sim_sq <- MicroSim(init_ppl, params, output, timesteps, agent_states, c_disc, PT.out = FALSE, strategy = "SQ", seed = calib.seed) # run for status quo
  return(c(
    colSums(matrix(sim_sq$oddeath, nrow = 12, ncol = 5))[-5],
    (colSums(matrix(sim_sq$fx_deaths, nrow = 12, ncol = 5)) / colSums(matrix(sim_sq$oddeath, nrow = 12, ncol = 5)))[-5],
    colSums(matrix(sim_sq$edvisits, nrow = 12, ncol = 5))[-5]
  ))
}
