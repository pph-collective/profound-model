#!/usr/bin/env Rscript

#' Main model entry point for the PROFOUND naloxone distribution model
#'
#' @description
#' `main()` sets up necessary parameters, runs the model, and saves the output
#'
#' @param TODO should change this to pretty much all be in main
#'
#' @returns
#' writes overdose stats to file
#'

# Load packages and scripts ---------------------------------------------------

rm(list = ls())
library("argparser")
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(tictoc)
source("parse_params.R")
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")

# just fixing some not in scope errors
if (FALSE) {
  parse_inputs <- NULL
  data_input <- NULL
  input_setup <- NULL
  output_setup <- NULL
  initiate_ppl <- NULL
  MicroSim <- NULL
}

args <- argparser::arg_parser("arguments")
args <- arg_parser::add_argument(args,
  "--inputs",
  help = "file containing input values",
  default = "params/base_params.yml"
)
argv <- argparser::parse_args(args)

inputs <- parse_inputs(argv$inputs)

# Set up data and param inputs/outputs
data <- data_input(inputs)
params <- input_setup(inputs, data)
params$data <- data
output <- output_setup(params)

if (file.exists(inputs$init_ppl_file)) { # import pop if possible
  init_ppl <- readRDS(inputs$init_ppl_file)
  print(paste0("Population loaded from file: ", inputs$init_ppl_file))
} else {
  init_ppl <- initiate_ppl(data, params$agent_states, seed = inputs$seed)
  saveRDS(init_ppl, inputs$init_ppl_file)
  print(paste0("Population saved to file: ", init_ppl.file))
}

main <- function(init_ppl, params, output) {
  print("start simulation")
  results_list <- data.frame()
  for (scenario in names(params$scenarios)) {
    results_list[scenario] <- MicroSim(init_ppl, params, output, scenario)
  }
  results <- do.call("rbind", results_list)
  write.csv(
    results,
    paste0(params$outdir, params$outfile),
    row.names = FALSE, na = "0"
  )
}

results <- data.frame(matrix(nrow = num_regions * 2, ncol = 6))
colnames(results) <- c("location", "scenario", "nlx_avail_rate", "nlx_avail", "overdose_deaths_rate", "overdose_deaths")

results$location <- rep(v.region, 2)
results$scenario <- rep(c("Status Quo", "Double"), each = length(v.region))
ppl_region <- colSums(Demographic[, -c(1:3)])

results$nlx_avail_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$avail_nlx) / ppl_region * 100000
results$nlx_avail[results$scenario == "Status Quo"] <- colSums(sim_sq$avail_nlx)
results$overdose_deaths_rate[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ]) / ppl_region * 100000
results$overdose_deaths[results$scenario == "Status Quo"] <- colSums(sim_sq$m.oddeath[49:60, ])
results$nlx_avail_rate[results$scenario == "Double"] <- colSums(sim_ep$avail_nlx) / ppl_region * 100000
results$nlx_avail[results$scenario == "Double"] <- colSums(sim_ep$avail_nlx)
results$overdose_deaths_rate[results$scenario == "Double"] <- colSums(sim_ep$m.oddeath[49:60, ] * 0.8) / ppl_region * 100000
results$overdose_deaths[results$scenario == "Double"] <- round(colSums(sim_ep$m.oddeath[49:60, ] * 0.8), 0)
write.csv(results, file = ("overdose_deaths.csv"))
