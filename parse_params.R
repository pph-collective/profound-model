#' parse inputs into a usable format
#'
#' This function takes in a yaml file of model options and converts it to a
#' named list. Params default to those stored in the `params/base_params.yml`
#' file.
#'
#' @param input_file The file containing the params.
#'
#' @return a named list with the script options
#'

library(yaml)
parse_inputs <- function(input_file) {
  if (FALSE) {
    get_item <- NULL
  }
  # read inputs and defaults
  inputs <- yaml::read_yaml(input_file)
  defaults <- yaml::read_yaml("params/base_params.yml")

  # create input list and add inputs to it
  input_params <- list()
  # general params
  input_params$year_start <- get_item(
    inputs$years$year_start,
    defaults$years$year_start
  )
  input_params$year_end <- get_item(
    inputs$years$year_end,
    defaults$years$year_end
  )
  input_params$discount <- get_item(inputs$discount$val, defaults$discount$val)
  input_params$seed <- get_item(inputs$seed$val, defaults$seed$val)
  input_params$strat <- get_item(inputs$strategy$val, defaults$strategy$val)
  input_params$scenarios <- get_item(inputs$scenarios, defaults$scenarios)

  # file and data locations
  input_params$main_table <- get_item(
    inputs$main_table$val,
    defaults$main_table$val
  )
  input_params$outdir <- get_item(inputs$out$dir, defaults$out$dir)
  input_params$outfile <- get_item(inputs$out$file, defaults$out$file)
  input_params$init_ppl_file <- get_item(
    inputs$init_ppl$val,
    defaults$init_ppl$val
  )

  # program data
  input_params$program_data <- get_item(
    inputs$program_data$val,
    defaults$program_data$val
  )

  # if random seed is desired, randomize
  if (input_params$seed == 0) {
    input_params$seed <- sample(1:2^15, 1)
  }

  # if the desired outfile doesn't exist, create it
  if (!dir.exists(input_params$outdir)) {
    dir.create(input_params$outdir, recursive = TRUE)
  }

  # write constructed params to file
  write_yaml(input_params, paste0(input_params$outdir, "params.yml"))
  return(input_params)
}


#' get an item or its default
#'
#' This function mimics the `get` option of python, but takes in an object
#' rather than object name, with a default if the object is null
#'
#' @param item The object to get.
#' @param alternate The default.
#'
#' @return The item or its alternate.
#'
get_item <- function(item, alternate) {
  if (is.null(item)) {
    return(alternate)
  }
  return(item)
}
