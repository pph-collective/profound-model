# TODO file header

parallel.fun <- function(calib.seed, params) {
  ## Run the simulation model over all sets of parameters
  sim_sq <- MicroSim(init_ppl, params, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = calib.seed) # run for status quo
  return(c(
    colSums(matrix(rowSums(sim_sq$m.oddeath.white), nrow = 12))[1:4],
    colSums(matrix(rowSums(sim_sq$m.oddeath.black), nrow = 12))[1:4],
    colSums(matrix(rowSums(sim_sq$m.oddeath.hisp), nrow = 12))[1:4],
    colSums(matrix(sim_sq$m.oddeath.fx.race[, "white"], nrow = 12))[1:4] / colSums(matrix(rowSums(sim_sq$m.oddeath.white), nrow = 12))[1:4],
    colSums(matrix(sim_sq$m.oddeath.fx.race[, "black"], nrow = 12))[1:4] / colSums(matrix(rowSums(sim_sq$m.oddeath.black), nrow = 12))[1:4],
    colSums(matrix(sim_sq$m.oddeath.fx.race[, "hisp"], nrow = 12))[1:4] / colSums(matrix(rowSums(sim_sq$m.oddeath.hisp), nrow = 12))[1:4],
    colSums(matrix(sim_sq$m.EDvisits, nrow = 12))
  ))
}
