#' Calculate costs at each timestep
#'
#' @description
#' `costs()` finds the total cost of interventions and the naloxone-specific
#' cost of the scenario.
#'
#' @param state Vector of current agent states.
#' @param ou_state Vector of current opioid use states.
#' @param nlx Total naloxone distributed at the timestep.
#' @param count Number of EMS and hospital visits in the timestep.
#' @param params Model parameters.
#'
#' @returns
#' A list with total cost at the timestep and naloxone-specific cost at the
#' timestep

costs <- function(state, ou_state, nlx, count_ems, count_hosp, params) {

  if (is.null(count_ems) && is.null(count_hosp)) {
    # if count is null, start from 0
    count_ems <- 0
    count_hosp <- 0
  }

  c.TC <- sum(state == "rx") * params$c_rx +
    sum(state == "il_lr") * params$c_il_lr +
    sum(state == "il_hr") * params$c_il_hr +
    sum(state == "inact") * params$c_inact +
    sum(state == "NODU") * params$c.NODU +
    sum(state == "relap" & ou_state == "rx") * params$c_relap["rx"] +
    sum(state == "relap" & ou_state == "il_lr") * params$c_relap["il_lr"] +
    sum(state == "relap" & ou_state == "il_hr") * params$c_relap["il_hr"] +
    count_ems * params$ems_cost +
    count_hosp * params$hosp_cost +
    nlx * (params$c_nlx_dtb + params$c_nlx_kit)

  c.nlx <- nlx * (params$c.nlx.dtb + params$c_nlx_kit)
  return(list(total_cost = c.TC, nlx_cost = c.nlx)) # return the costs
}
