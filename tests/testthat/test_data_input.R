library(testthat)

source("data_input.R")

#' Test the data input function
#'
#' @description
#' Tests the capabilities of the data input functions
#'
#' * `test_data_input()` tests that all parameters are as expected
#'
#'

test_data_input <- function() {
  params <- parse_inputs("tests/params/test_params.yml")
  data <- data_input(params$main_table)

  # initial demographics
  expect_equal(1000, data$initials$ppl_size)
  expect_equal(0.68, data$initials$prev.oud)
  expect_equal(0.32, data$initials$prev.NODU.m)
  expect_equal(0.32, data$initials$prev.NODU.f)

  # race, age, region
  expect_length(data$v.region, 39)
  age <- rep(
    c("12to17", "18to25", "26to34", "35to49", "50to64", "65older"),
    8
  )
  expect_identical(data$v.demo.age, age)

  # check the correct races are present
  races <- count(as.data.frame(data$v.demo.race), data$v.demo.race)
  expect_length(races$n, 4)
  for (i in races[, 2]) {
    expect_equal(i, 12)
  }

  # note: equations from excel doc main_table
  # illicit
  expect_equal(data$initials$ini.il.m, (93000+62000)/(93000+62000+781000))
  expect_equal(data$initials$ini.il.f, (30000+55000)/(30000+55000+616000))

  # high risk among illicit
  expect_equal(.486, data$initials$ini.il.hr.m)
  expect_equal(data$initials$ini.il.hr.f, 0.4105)

  # fentanyl exposure
  expect_equal(data$init_oud_fx, 22/39)

  # stimulant
  expect_equal(data$ini.NOUD.fx, 0.1567, tolerance = 0.00001)
  expect_equal(data$initials$ini.everod.sti, 0.252445, tolerance = 0.00001)

  # mortality
  mort <- 1 - (1 - 13.3 / 1000000)^ (1/12)
  expect_equal(data$mor.bg[1], 1 - (1 - 13.3 / 1000000)^ (1/12))

  # overdose risk
  expect_equal(data$od.preb.sub, 0.007543849358)

  # transition
  expect_equal(data$p.preb2il.lr, 0.000418)

  # overdose ems
  expect_equal(data$OD_loc[1], 0.692, tolerance = 0.001)

  # decision tree
  expect_equal(data$OD_wit_priv, 0.645)

  # mortality
  expect_equal(data$mor_bl, 0.1111111)

  # naloxone
  woonsocket <- data.frame(high = c(34, 112, 430, 584, 667))
  woonsocket$low <- c(0, 0, 1, 71, 170)
  rownames(woonsocket) <- 2015:2019
  rownames(woonsocket) <- as.character(rownames(woonsocket))
  x <- as.data.frame(data$NxOEND[, , "Woonsocket"])
  expect_equal(x, woonsocket)

  # Cost
  expect_equal(data$c.preb, 1019)
}
