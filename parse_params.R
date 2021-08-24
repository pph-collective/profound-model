#'
#' 
library(yaml)


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
parse_inputs <- function(input_file){
    # read in inputs and default
    inputs <- read_yaml(input_file)
    defaults <- read_yaml("params/base_params.yml")

    # create input list and add inputs to it
    input_params <- list()
    input_params$year_start <- get_item(inputs$years$year_start, defaults$years$year_start)
    input_params$year_end <- get_item(inputs$years$year_end, defaults$years$year_end)
    input_params$main_table <- get_item(inputs$main_table$val, defaults$main_table$val)
    input_params$discount <- get_item(inputs$discount$val, defaults$input$val)
    input_params$seed <- get_item(inputs$seed$val, defaults$seed$val)
    input_params$strat <- get_item(inputs$strategy$val, defaults$strategy$val)
    input_params$outdir <- get_item(inputs$outfolder$val, defaults$outfile$val)
    # put the output file in the outdir
    input_params$outfile <- paste0(input_params$outdir, "overdoses.csv")
    input_params$init_ppl_file <- get_item(inputs$init_ppl$val, inputs$init_ppl$val)
    input_params$scenarios <- inputs$scenarios

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
    if (is.null(item)) { return(alternate) }
    return(item)
}

