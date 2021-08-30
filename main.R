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

# Load packages and scripts -----------------------------------------------------

library("argparser")
library(dplyr)
library(tictoc)
library(openxlsx)
library(abind)
library(tictoc)
source("population_initialization.R")
source("transition_probability.R")
source("microsim.R")
source("decision_tree.R")
source("data_input.R")
source("naloxone_availability.R")
source("cost_effectiveness.R")
source("parse_params.R")
source("io_setup.R")
source("data_input.R")

# parse command line arguments--------------------------------------------------
args <- arg_parser("arguments")
args <- add_argument(args, "--inputs", help = "file containing input values", default = "params/base_params.yml")
argv <- parse_args(args)

inputs <- parse_inputs(argv$inputs)
d.c <- inputs$discount # discounting of costs by 3%
main_table <- inputs$main_table

sw.EMS.ODloc <- inputs$strat
out.file <- inputs$outfile
# init_ppl.file <- inputs$init_ppl


data <- data_input(main_table)  # empirical
# add parameters
params <- input_setup(inputs, data)
# create output table
output <- output_setup(params)

## Initialize the study population - people who are at risk of opioid overdose
ppl_info <- c(
  "sex",
  "race",
  "age",
  "residence",
  "curr.state",
  "OU.state",
  "init.age",
  "init.state",
  "ever.od",
  "fx"
)

if (file.exists(inputs$init_ppl_file)) { # import pop if possible
  init_ppl <- readRDS(inputs$init_ppl_file)
  print(paste0("Population loaded from file: ", inputs$init_ppl_file))
} else { # otherwise, create pop
  init_ppl <- initiate_ppl(data, seed = inputs$seed)
  saveRDS(init_ppl, inputs$init_ppl_file)
  print(paste0("Population saved to file: ", inputs$init_ppl_file))
}


# Run the simulation =============================

main <- function(init_ppl, params, data, output, timesteps, d.c, expansion, seed){
  # what I want to do: for each scenario, run the simulation. If it's a program, send it to the program eval to change the probs. Otherwise, scale as needed
  print("simulate")
  results <- data.frame()
  for (scenario in names(params$scenarios)){
    results <- rbind(results, MicroSim(init_ppl, params, data, output, d.c, scenario, seed))
    # outfile <- paste0(params$outdir, scenario, "_overdose.csv")
    outfile <- "overdoses.csv"
    write.csv(results, paste0(params$outdir, params$outfile), row.names = FALSE, na="0")
  }
}

main(init_ppl, params, data, output, timesteps, d.c, expansion, inputs$seed)

