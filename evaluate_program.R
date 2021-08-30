###########################################################################################################
#####################            Program evaluation - individual program             ######################
###########################################################################################################
# To evaluate the impact of specific OEND program at a few hypothetical expansion levels

#############################################################################
# 1. SET directory and workspace
#############################################################################
# library(openxlsx)
# library(dplyr)
# library(abind)
# source("population_initialization.R")
# source("transition_probability.R")
# source("microsim.R")
# source("decision_tree.R")
# source("data_input.R")
# source("naloxone_availability.R")
# source("cost_effectiveness.R")


evaluate_program <- function(params) {# Program data
  program_data <- read.xlsx(params$program_data$program_file$val, sheet = params$program_data$val)
  program_levels <- c(1, 5, 10, 20, 50)
  program_add.array <- array(0, dim = c(length(program_levels), 2, dim(program_data)[1]))
  # REVIEWED add.array = additional naloxone for high or low risk program
  dimnames(program_add.array)[[2]] <- c("high", "low")
  for (i in 1:length(program_levels)) {
    if (program_data$Risk[1] == "high") {
      program_add.array[i, "high", ] <- round(program_data$Volume[1] * program_data$Proportion * program_levels[i], 0)
      program_add.array[i, "low", ] <- 0
    } else {
      program_add.array[i, "high", ] <- 0
      program_add.array[i, "low", ] <- round(program_data$Volume[1] * program_data$Proportion * program_levels[i], 0)
    }
  }

  # ## Initialize the study population - people who are at risk of opioid overdose
  # ppl_info <- c("sex", "race", "age", "residence", "curr.state", "OU.state", "init.age", "init.state", "ever.od", "fx")
  # if (file.exists(paste0("Inputs/init_pop.rds"))) {
  #   init_ppl <- readRDS(paste0("Inputs/init_pop.rds"))
  # } else if (!file.exists(paste0("Inputs/init_pop.rds"))) {
  #   init_ppl <- ppl.initiation(initials = initials, seed = seed)
  #   saveRDS(init_ppl, paste0("Inputs/init_pop.rds"))
  # }

  # REVIEWED ls = list
  simulation_data <- readRDS(file = paste0("calibration/CalibratedData.rds"))
  simulation_seed <- readRDS(file = paste0("calibration/CalibratedSeed.rds"))
  simulation_seed <- simulation_seed[1:100]

  # REVIEWED sq = status quo; pg = program
  sq.dh.mx <- sq.nx.mx <- matrix(0, nrow = length(v.region), ncol = length(simulation_seed))
  program_dh.ar <- program_nx.ar <- array(0, dim = c(dim(program_add.array)[1], length(v.region), length(simulation_seed)))
  nlx.used.mx <- matrix(0, nrow = length(simulation_seed), ncol = 1 + length(program_levels))
  od.death.mx <- matrix(0, nrow = length(simulation_seed), ncol = 1 + length(program_levels))
  scenario.name <- c("Status Quo", "100% increase", "500% increase", "1000% increase", "2000% increase", "5000% increase")
  colnames(nlx.used.mx) <- scenario.name
  colnames(od.death.mx) <- scenario.name

  for (ss in 1:length(simulation_seed)) {
    print(paste0("Parameter set: ", ss))
    params.temp <- simulation_data[[ss]]
    params.temp$NxDataPharm$pe <- 0
    params.temp$mortality_nx <- params.temp$mor_bl * (1 - 0.9)
    sim_sq <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, inputs$discount, PT.out = FALSE, strategy = "SQ", seed = simulation_seed[ss]) # run for status quo
    sq.dh.mx[, ss] <- colSums(sim_sq$m.oddeath[(timesteps - 11):timesteps, ])
    sq.nx.mx[, ss] <- colSums(sim_sq$n.nlx.OEND.str)
    nlx.used.mx[ss, "Status Quo"] <- sum(sim_sq$v.nlxused[(timesteps - 11):timesteps])
    od.death.mx[ss, "Status Quo"] <- sum(sim_sq$m.oddeath[(timesteps - 11):timesteps, ])

    for (ll in 1:dim(program_add.array)[1]) {
      params.temp$program_add <- program_add.array[ll, , ]
      sim_pg <- MicroSim(init_ppl, params = params.temp, timesteps, agent_states, inputs$discount, PT.out = FALSE, strategy = "program", seed = simulation_seed[ss]) # run for program scenario
      program_dh.ar[ll, , ss] <- colSums(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])
      program_nx.ar[ll, , ss] <- colSums(sim_pg$n.nlx.OEND.str)
      nlx.used.mx[ss, scenario.name[ll + 1]] <- sum(sim_pg$v.nlxused[(timesteps - 11):timesteps])
      od.death.mx[ss, scenario.name[ll + 1]] <- sum(sim_pg$m.oddeath[(timesteps - 11):timesteps, ])
    }
  }

  pop.rgn <- colSums(params$demographic[, -c(1:3)]) # population in each region (for calculation of rate)

  NoDeaths <- data.frame(matrix(nrow = n.rgn * (1 + dim(program_add.array)[1]), ncol = 5))
  x <- c("location", "scenario", "mean", "upper", "lower")
  colnames(NoDeaths) <- x
  NoDeaths$location <- rep(v.rgn, 1 + dim(program_add.array)[1])

  NoDeaths$scenario <- rep(scenario.name, each = length(v.rgn))

  RateNlx <- NoNlx <- RateDeaths <- NoDeaths
  # Number of deaths
  NoDeaths$mean[NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, mean)
  NoDeaths$upper[NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, quantile, probs = 0.975)
  NoDeaths$lower[NoDeaths$scenario == "Status Quo"] <- apply(sq.dh.mx, 1, quantile, probs = 0.025)

  for (sc in 2:length(scenario.name)) {
    NoDeaths$mean[NoDeaths$scenario == scenario.name[sc]] <- apply(program_dh.ar[sc - 1, , ], 1, mean)
    NoDeaths$upper[NoDeaths$scenario == scenario.name[sc]] <- apply(program_dh.ar[sc - 1, , ], 1, quantile, probs = 0.975)
    NoDeaths$lower[NoDeaths$scenario == scenario.name[sc]] <- apply(program_dh.ar[sc - 1, , ], 1, quantile, probs = 0.025)
  }

  # Rate of deaths
  RateDeaths[, c("mean", "upper", "lower")] <- NoDeaths[, c("mean", "upper", "lower")] / pop.rgn * 100000

  # Number of Naloxone kits
  NoNlx$mean[NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, mean)
  NoNlx$upper[NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, quantile, probs = 0.975)
  NoNlx$lower[NoNlx$scenario == "Status Quo"] <- apply(sq.nx.mx, 1, quantile, probs = 0.025)

  for (sc in 2:length(scenario.name)) {
    NoNlx$mean[NoNlx$scenario == scenario.name[sc]] <- apply(program_nx.ar[sc - 1, , ], 1, mean)
    NoNlx$upper[NoNlx$scenario == scenario.name[sc]] <- apply(program_nx.ar[sc - 1, , ], 1, quantile, probs = 0.975)
    NoNlx$lower[NoNlx$scenario == scenario.name[sc]] <- apply(program_nx.ar[sc - 1, , ], 1, quantile, probs = 0.025)
  }

  # Rate of Naloxone kits
  RateNlx[, c("mean", "upper", "lower")] <- NoNlx[, c("mean", "upper", "lower")] / pop.rgn * 100000

  write.csv(NoDeaths, file = ("Outputs/Program/Number.Deaths.csv"), row.names = F)
  write.csv(RateDeaths, file = ("Outputs/Program//Rate.Deaths.csv"), row.names = F)
  write.csv(NoNlx, file = ("Outputs/Program//Number.Naloxone.csv"), row.names = F)
  write.csv(RateNlx, file = ("Outputs/Program//Rate.Naloxone.csv"), row.names = F)
  write.csv(nlx.used.mx, file = ("Outputs/Program//NaloxoneUseinputs$discountsv"), row.names = F)
  write.csv(od.death.mx, file = ("Outputs/Program//TotalODdeaths.csv"), row.names = F)
}
