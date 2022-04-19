# TODO file header

parallel.fun <- function(calib.seed, params) {
  ## Run the simulation model over all sets of parameters
  sim_sq <- MicroSim(init_ppl, params, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = calib.seed) # run for status quo
  return(c(
    colSums(matrix(sim_sq$v.oddeath, nrow = 12)),
    (colSums(matrix(sim_sq$m.oddeath.fx, nrow = 12)) / colSums(matrix(sim_sq$v.oddeath, nrow = 12))),
    colSums(matrix(sim_sq$m.EDvisits, nrow = 12))
  ))
}
