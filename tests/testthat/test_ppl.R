library(testthat)

source("population_initialization.R")
source("data_input.R")

test_init_ppl <- function() {
    main_table <- "tests/params/main_table.xlsx"
    states <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")
    data <- data_input(main_table)
    ppl <- initiate_ppl(data, states)
    # import workbook for comparison
    WB <- loadWorkbook(main_table)
    init_ppl_params <- read.xlsx(WB, sheet = "InitialPop")
    expect_equal(nrow(ppl), init_ppl_params$pe[init_ppl_params$par == "ppl_size"])
    print("tests pass!")
}
