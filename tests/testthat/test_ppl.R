library(testthat)

source("population_initialization.R")
source("data_input.R")

test_init_ppl <- function() {
    main_table <- "tests/params/main_table.xlsx"
    states <- c("rx", "il_lr", "il_hr", "inact", "NODU", "relap", "dead")
    data <- data_input(main_table)
    ppl <- initiate_ppl(data, states)
    # import workbook for comparison
    workbook <- loadWorkbook(main_table)
    init_ppl_params <- read.xlsx(workbook, sheet = "InitialPop")
    expect_equal(nrow(ppl), init_ppl_params$pe[init_ppl_params$par == "ppl_size"])
    print("tests pass!")
}
