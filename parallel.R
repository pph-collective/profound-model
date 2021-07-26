# TODO file header

parallel.fun <- function(calib.seed, params) {
  ## Run the simulation model over all sets of parameters
  sim_sq <- MicroSim(init_ppl, params, timesteps, agent_states, d.c, PT.out = FALSE, strategy = "SQ", seed = calib.seed) # run for status quo
  return(c(
    colSums(matrix(sim_sq$v.oddeath, nrow = 12, ncol = 5))[-5],
    (colSums(matrix(sim_sq$m.oddeath.fx, nrow = 12, ncol = 5)) / colSums(matrix(sim_sq$v.oddeath, nrow = 12, ncol = 5)))[-5],
    colSums(matrix(sim_sq$m.EDvisits, nrow = 12, ncol = 5))[-5]
  ))
}
