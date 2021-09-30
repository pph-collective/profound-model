#' Calculate costs at each timestep
#'
#' @description
#' `Costs()` finds the total cost of interventions and the naloxone-specific costs
#' of the scenario.
#'
#' @param state Vector of current agent states.
#' @param OU.state Vector of current opioid use states.
#' @param nlx Total naloxone distributed at the timestep.
#' @param count Number of EMS and hospital visits in the timestep.
#' @param params Model parameters.
#'
#' @returns
#' A list with total cost at the timestep and naloxone-specific cost at the timestep

Costs <- function(state, OU.state, nlx, count, params) {

  if (is.null(count)) {
    count.EMS <- 0
    count.hospcare <- 0
  } else {
    count.EMS <- params$count$n.EMS
    count.hospcare <- params$count$n.hospcare
  }
  c.TC <- sum(state == "preb") * params$c.preb +
    sum(state == "il.lr") * params$c.il.lr +
    sum(state == "il.hr") * params$c.il.hr +
    sum(state == "inact") * params$c.inact +
    sum(state == "NODU") * params$c.NODU +
    sum(state == "relap" & OU.state == "preb") * params$c.relap.v["preb"] +
    sum(state == "relap" & OU.state == "il.lr") * params$c.relap.v["il.lr"] +
    sum(state == "relap" & OU.state == "il.hr") * params$c.relap.v["il.hr"] +
    count.EMS * params$c.EMS +
    count.hospcare * params$c.hospcare +
    nlx * (params$c.nlx.dtb + params$c.nlx.kit)

  c.nlx <- nlx * (params$c.nlx.dtb + params$c.nlx.kit)
  return(list(total_cost = c.TC, nlx_cost = c.nlx)) # return the costs
}
