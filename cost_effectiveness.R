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
  list2env(params, environment())
  if (is.null(count)) {
    count.EMS <- 0
    count.hospcare <- 0
  } else {
    count.EMS <- count$n.EMS
    count.hospcare <- count$n.hospcare
  }
  c.TC <- sum(state == "preb") * c.preb +
    sum(state == "il.lr") * c.il.lr +
    sum(state == "il.hr") * c.il.hr +
    sum(state == "inact") * c.inact +
    sum(state == "NODU") * c.NODU +
    sum(state == "relap" & OU.state == "preb") * c.relap.v["preb"] +
    sum(state == "relap" & OU.state == "il.lr") * c.relap.v["il.lr"] +
    sum(state == "relap" & OU.state == "il.hr") * c.relap.v["il.hr"] +
    count.EMS * c.EMS +
    count.hospcare * c.hospcare +
    nlx * (c.nlx.dtb + c.nlx.kit)

  c.nlx <- nlx * (c.nlx.dtb + c.nlx.kit)
  return(list(total_cost = c.TC, nlx_cost = c.nlx)) # return the costs
}
