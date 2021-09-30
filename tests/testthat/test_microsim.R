library(testthat)
library(abind)
source("io_setup.R")
source("microsim.R")
source("parse_params.R")
source("io_setup.R")
source("data_input.R")

test_step_one <- function() {
    init_ppl <- readRDS("tests/init_ppl.RDS")
    params <- parse_inputs("tests/params/test_params_empty.yml")
    data <- data_input(params$main_table)
    params <- input_setup(params, data)
    params$timesteps <- 1
    output <- output_setup(params)
    ppl <- list()
    ppl[[1]] <- init_ppl
    
    # do not include pharmacy data
    nlx_array <- data$NxOEND
    nlx_array <- abind(nlx_array, nlx_array[dim(nlx_array)[1], , ], along = 1)
    x <- step(1, output, nlx_array, ppl, data, 1, params, 0)
    ppl_2 <- x$ppl_list[[1]]
    ppl_1 <- ppl[[1]]

    # full list should be the same except fentanyl
    expect_identical(ppl_1[1, -which(names(ppl_1) == "fx")], ppl_2[1, -which(names(ppl_2) == "fx")])
}
