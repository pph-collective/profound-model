################################################################################################
#########################       Function for parallel calibration     ##########################
################################################################################################

################################################################################################
# Put simulation function and returning required outcomes in one function
################################################################################################

parallel.fun <- function(calib.seed, vparameters){
  ## Run the simulation model over all sets of parameters
  sim_sq    <- MicroSim(init.pop, vparameters, n.t, v.state, d.c, PT.out = FALSE, Str = "SQ", seed = calib.seed)        # run for status quo
  return(c(colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5))[-5], 
           (colSums(matrix(sim_sq$m.oddeath.fx, nrow=12, ncol=5)) / colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5)))[-5],
           colSums(matrix(sim_sq$m.EDvisits, nrow=12, ncol=5))[-5]))
}