###############################################################################################
#######################         City-level validation         #################################
###############################################################################################
# Module for comparing model projections at the city level with observed surveillance data
#
# Authors: Xiao Zang, PhD, Sam Bessey, MS
#
# People, Place and Health Collective, Department of Epidemiology, Brown University
#
###############################################################################################

#############################################################################
# 1. SET directory and workspace
#############################################################################
rm(list = ls())
library(argparser)
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(tictoc)
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")

yr_start <- 2016 # starting year of simulation
yr_end <- 2020 # end year of simulation (also the year for evaluation)
discount.rate <- 0.03 # discounting of costs by 3%

args <- arg_parser("arguments")
args <- add_argument(args, "--seed", help = "seed for random numbers", default = 2021)
args <- add_argument(args, "--regional", help = "flag to run regional model", flag = TRUE)
args <- add_argument(args, "--outfile", help = "file to store outputs", default = "OverdoseDeath_RIV1_0.csv")
args <- add_argument(args, "--ppl", help = "file with initial ppl info", default = "Inputs/init_pop.rds")
argv <- parse_args(args)
seed <- as.integer(argv$seed)
init_ppl.file <- argv$ppl
source("io_setup.R")

##################################### Run simulation with calibrated parameters ######################################################
sim.seed <- sim.seed[1:500]
ODdeaths16 <- matrix(0, nrow = num_regions, ncol = length(sim.seed))
ODdeaths17 <- matrix(0, nrow = num_regions, ncol = length(sim.seed))
ODdeaths18 <- matrix(0, nrow = num_regions, ncol = length(sim.seed))
ODdeaths19 <- matrix(0, nrow = num_regions, ncol = length(sim.seed))
row.names(ODdeaths16) <- row.names(ODdeaths17) <- row.names(ODdeaths18) <- row.names(ODdeaths19) <- v.region
for (ss in 1:length(sim.seed)) {
  print(paste0("Parameter set: ", ss))
  vparameters.temp <- sim.data.ls[[ss]]
  sim_sq <- MicroSim(init_ppl, params = vparameters.temp, timesteps, agent_states, discount.rate, PT.out = FALSE, strategy = "SQ", seed = sim.seed[ss])        # run for status quo
  ODdeaths16[, ss] <- colSums(sim_sq$m.oddeath[1:12, ])
  ODdeaths17[, ss] <- colSums(sim_sq$m.oddeath[13:24, ])
  ODdeaths18[, ss] <- colSums(sim_sq$m.oddeath[25:36, ])
  ODdeaths19[, ss] <- colSums(sim_sq$m.oddeath[37:48, ])
}

detach("package:openxlsx", unload = TRUE)
library(xlsx)

write.xlsx(ODdeaths16,
  file = "CityLevelValidation.xlsx", sheetName = "2016",
  col.names = F, row.names = T
)
write.xlsx(ODdeaths17,
  file = "CityLevelValidation.xlsx", sheetName = "2017", append = TRUE,
  col.names = F, row.names = T
)
write.xlsx(ODdeaths18,
  file = "CityLevelValidation.xlsx", sheetName = "2018", append = TRUE,
  col.names = F, row.names = T
)
write.xlsx(ODdeaths19,
  file = "CityLevelValidation.xlsx", sheetName = "2019", append = TRUE,
  col.names = F, row.names = T
)


