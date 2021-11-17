library(testthat)

source("parse_params.R")
source("population_initialization.R")
source("data_input.R")

test_init_ppl <- function() {
  # main_table <- "tests/params/main_table.xlsx"
  states <- c("rx", "il_lr", "il_hr", "inact", "NODU", "relap", "dead")
  params <- parse_inputs("tests/params/test_params.yml")
  data <- data_input(params)
  environment(initiate_ppl) <- environment(data_input)
  ppl <- initiate_ppl(data, states)
  # import workbook for comparison
  workbook <- loadWorkbook(params$main_table)
  init_ppl_params <- read.xlsx(workbook, sheet = "InitialPop")
  expect_equal(
    nrow(ppl),
    init_ppl_params$pe[init_ppl_params$par == "pop.size"]
  )
  print("tests pass!")
}

