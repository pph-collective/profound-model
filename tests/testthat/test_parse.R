library(testthat)
library(yaml)

source("parse_params.R")


test_get_item <- function() {
    expect_identical(1, get_item(NULL, 1))
    expect_identical(5, get_item(5, 1))
    expect_identical(0, get_item(0, 1))
}

test_parse_empty <- function() {
    base <- parse_inputs("params/base_params.yml")
    params <- parse_inputs("tests/params/test_params_empty.yml")
    multi <- parse_inputs("tests/params/test_multi_scenario.yml")

    expect_identical(params, base)
    # everything but scenarios should be the same
    for (param in names(params)) {
        if (param != "scenarios") {
            expect_identical(params[param], multi[param])
        }
    }
    # multi scenario has two scenarios
    expect_equal(length(multi$scenarios), 2)
    expect_identical(multi$scenarios$SQ, params$scenario$SQ)
}

test_new_scenario <- function() {
    expand <- parse_inputs("tests/params/test_expansion.yml")
    # should not have status quo scenario
    expect_length(expand$scenarios, 1)
    expect_equal(names(expand$scenarios), "expand")
}
