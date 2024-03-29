###########################################################################################
#########################       Input and Output setup     ################################
###########################################################################################

###########################################################################################
# This is to define all the input and output variables required in the model (consistent across all scenarios)

# INPUT setup
pop.info <- c(
  "sex", "race", "age", "residence", "curr.state",
  "OU.state", "init.age", "init.state", "ever.od", "fx"
) # information for each model individual
agent_states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead") # vector for state names
v.oustate <- c("preb", "il.lr", "il.hr") # vector for active opioid use state names
num_states <- length(agent_states) # number of states
num_years <- yr_end - yr_start + 1
timesteps <- 12 * num_years # number of time cycles (in month)
num_regions <- length(v.region) # number of regions

# Load or generate initial study population: people who are at risk of opioid overdose
ppl_info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
if (file.exists(init_ppl.file)) {
  init_ppl <- readRDS(init_ppl.file)
  print(paste0("Population loaded from file: ", init_ppl.file))
} else {
  source("population_initialization.R")
  init_ppl <- initiate_ppl(initials = initials, seed = seed)
  saveRDS(init_ppl, init_ppl.file)
  print(paste0("Population saved to file: ", init_ppl.file))
}

# OUTPUT matrices and vectors
v.od <- rep(0, times = timesteps) # count of overdose events at each time step
v.oddeath <- rep(0, times = timesteps) # count of overdose deaths at each time step
v.oddeath.w <- rep(0, times = timesteps) # count of overdose deaths that were witnessed at each time step
m.oddeath <- matrix(0, nrow = timesteps, ncol = num_regions)
colnames(m.oddeath) <- v.region
v.odpriv <- rep(0, times = timesteps) # count of overdose events occurred at private setting at each time step
v.odpubl <- rep(0, times = timesteps) # count of overdose events occurred at public setting at each time step
v.deathpriv <- rep(0, times = timesteps) # count of overdose deaths occurred at private setting at each time step
v.deathpubl <- rep(0, times = timesteps) # count of overdose deaths occurred at public setting at each time step
v.nlxused <- rep(0, times = timesteps) # count of naloxone kits used at each time step
v.str <- c("SQ", "expand", "program") # store the strategy names
cost.item <- c("TotalCost", "NxCost")
cost.matrix <- matrix(0, nrow = timesteps, ncol = length(cost.item))
colnames(cost.matrix) <- cost.item
m.oddeath.fx <- rep(0, times = timesteps) # count of overdose deaths with fentanyl present at each time step
m.oddeath.op <- rep(0, times = timesteps) # count of overdose deaths among opioid users at each time step
m.oddeath.st <- rep(0, times = timesteps) # count of overdose deaths among stimulant users at each time step
m.EDvisits <- rep(0, times = timesteps) # count of opioid overdose-related ED visits at each time step
m.oddeath.hr <- rep(0, times = timesteps) # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
m.oddeath.preb <- m.oddeath.il.lr <- m.oddeath.il.hr <- m.nlx.mn.OEND <- m.nlx.mn.Pharm <- rep(0, times = timesteps)   # count of overdose deaths stratified by risk group
# v.nlx.mn <- matrix(0, nrow= length(sim.seed), ncol = timesteps)  # count of naloxone kits in circulation in each month

## Load calibrated parameters and seeds
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))