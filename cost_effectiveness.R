########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Module for Data Input of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# People, Place and Health Collective, Department of Epidemiology, Brown University
#
# Created: Dec 30, 2020
# Last update: Jan 31, 2021
#
##################################################################################
###################        Cost effectiveness analysis      ######################
##################################################################################

### Costs function
# The Costs function estimates the costs at every cycle.

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
