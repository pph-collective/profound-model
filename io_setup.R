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
  ## TEMPORARY: remove "other" population 
  init_ppl <- subset(init_ppl, race != "other")
  init_ppl$ind <- 1:nrow(init_ppl)
  saveRDS(init_ppl, init_ppl.file)
  print(paste0("Population saved to file: ", init_ppl.file))
}

# # (NEW) Characterize the entry of population into the model (fixed number each month and randomly drawn from the initial population)
# add_pop_file <- "Inputs/add_pop.rds"
# if(file.exists(add_pop_file)) {
#   add_pop <- readRDS(add_pop_file)
# } else {
#   set.seed(2023)
#   # init_HS <- init_ppl[init_ppl$OU.state == "il.hr" | init_ppl$OU.state == "il.lr" | init_ppl$OU.state == "NODU", ]
#   init_HS <- init_ppl[init_ppl$OU.state == "il.hr" | init_ppl$OU.state == "il.lr", ]
#   add_pop <- init_HS[sample(1:nrow(init_HS), annual.entry * (yr_end - yr_start + 1), replace = T), ]
#   add_pop$ind <- c((last(init_ppl$ind)+1):(last(init_ppl$ind)+nrow(add_pop)))
#   add_pop$ever.od <- 0
#   row.names(add_pop) <- NULL
#   saveRDS(add_pop, add_pop_file)
# }

# OUTPUT matrices and vectors
v.od <- rep(0, times = timesteps) # count of overdose events at each time step
v.oddeath <- rep(0, times = timesteps) # count of overdose deaths at each time step
v.oddeath.w <- rep(0, times = timesteps) # count of overdose deaths that were witnessed at each time step
v.oddeath.w.race <- matrix(0, nrow = timesteps, ncol = 3)
colnames(v.oddeath.w.race) <- c("white", "black", "hisp")
m.oddeath <- matrix(0, nrow = timesteps, ncol = num_regions)
m.oddeath.white <- matrix(0, nrow = timesteps, ncol = num_regions)
m.oddeath.black <- matrix(0, nrow = timesteps, ncol = num_regions)
m.oddeath.hisp  <- matrix(0, nrow = timesteps, ncol = num_regions)
colnames(m.oddeath) <- colnames(m.oddeath.white) <- colnames(m.oddeath.black) <- colnames(m.oddeath.hisp) <- v.region
v.odpriv <- rep(0, times = timesteps) # count of overdose events occurred at private setting at each time step
v.odpubl <- rep(0, times = timesteps) # count of overdose events occurred at public setting at each time step
v.deathpriv <- rep(0, times = timesteps) # count of overdose deaths occurred at private setting at each time step
v.deathpubl <- rep(0, times = timesteps) # count of overdose deaths occurred at public setting at each time step
v.nlxused <- rep(0, times = timesteps) # count of naloxone kits used at each time step
v.str <- c("SQ", "non-target", "equity") # store the strategy names
cost.item <- c("TotalCost", "NxCost")
cost.matrix <- matrix(0, nrow = 3*12, ncol = length(cost.item))
colnames(cost.matrix) <- cost.item
qaly.item <- c("all", "white", "black", "hisp")
qaly.matrix <- matrix(0, nrow = 3*12, ncol = length(qaly.item))
colnames(qaly.matrix) <- qaly.item
m.oddeath.fx <- rep(0, times = timesteps) # count of overdose deaths with fentanyl present at each time step
m.oddeath.op <- rep(0, times = timesteps) # count of overdose deaths among opioid users at each time step
m.oddeath.st <- rep(0, times = timesteps) # count of overdose deaths among stimulant users at each time step
m.EDvisits <- rep(0, times = timesteps) # count of opioid overdose-related ED visits at each time step
m.oddeath.hr <- rep(0, times = timesteps) # count of overdose deaths among high-risk opioid users (inject heroin) at each time step
m.oddeath.preb <- m.oddeath.il.lr <- m.oddeath.il.hr <- m.nlx.mn.OEND <- m.nlx.mn.Pharm <- rep(0, times = timesteps)   # count of overdose deaths stratified by risk group
m.oddeath.fx.race <- m.EDvisits.race <- matrix(0, nrow = timesteps, ncol = 3)
colnames(m.oddeath.fx.race) <- colnames(m.EDvisits.race) <- c("white", "black", "hisp")
# v.nlx.mn <- matrix(0, nrow= length(sim.seed), ncol = timesteps)  # count of naloxone kits in circulation in each month

## Load calibrated parameters and seeds
sim.data.ls <- readRDS(file = paste0("calibration/CalibratedData.rds"))
sim.seed    <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))