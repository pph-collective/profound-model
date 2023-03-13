########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Module for Data Input of the Profound Naloxone distribution model:
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# People, Place and Health Collective, Department of Epidemiology, Brown University
#
# Created: Dec 30, 2020
# Last update: Jan 25, 2022
#
##################################################################################
###################        Cost effectiveness analysis      ######################
##################################################################################

### Costs function
# The Costs function estimates the costs and qalys at every cycle.

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
  return(c(c.TC, c.nlx)) # return the costs
}

QALYs <- function(state, OU.state, race, params) {
  list2env(params, environment())
  # if (is.null(count)) {
  #   count.EMS <- 0
  #   count.hospcare <- 0
  # } else {
  #   count.EMS <- count$n.EMS
  #   count.hospcare <- count$n.hospcare
  # }
  Q.TQ <- 
    sum(state == "preb") * q.preb +
    sum(state == "il.lr") * q.il.lr +
    sum(state == "il.hr") * q.il.hr +
    sum(state == "inact") * q.inact +
    sum(state == "NODU") * q.NODU +
    sum(state == "relap" & OU.state == "preb") * q.relap.v["preb"] +
    sum(state == "relap" & OU.state == "il.lr") * q.relap.v["il.lr"] +
    sum(state == "relap" & OU.state == "il.hr") * q.relap.v["il.hr"] 
  
  Q.white <- 
    sum(state == "preb" & race == "white") * q.preb +
    sum(state == "il.lr" & race == "white") * q.il.lr +
    sum(state == "il.hr" & race == "white") * q.il.hr +
    sum(state == "inact" & race == "white") * q.inact +
    sum(state == "NODU" & race == "white") * q.NODU +
    sum(state == "relap" & OU.state == "preb" & race == "white") * q.relap.v["preb"] +
    sum(state == "relap" & OU.state == "il.lr" & race == "white") * q.relap.v["il.lr"] +
    sum(state == "relap" & OU.state == "il.hr" & race == "white") * q.relap.v["il.hr"]
  
  Q.black <- 
    sum(state == "preb" & race == "black") * q.preb +
    sum(state == "il.lr" & race == "black") * q.il.lr +
    sum(state == "il.hr" & race == "black") * q.il.hr +
    sum(state == "inact" & race == "black") * q.inact +
    sum(state == "NODU" & race == "black") * q.NODU +
    sum(state == "relap" & OU.state == "preb" & race == "black") * q.relap.v["preb"] +
    sum(state == "relap" & OU.state == "il.lr" & race == "black") * q.relap.v["il.lr"] +
    sum(state == "relap" & OU.state == "il.hr" & race == "black") * q.relap.v["il.hr"]
  
  Q.hisp <- 
    sum(state == "preb" & race == "hisp") * q.preb +
    sum(state == "il.lr" & race == "hisp") * q.il.lr +
    sum(state == "il.hr" & race == "hisp") * q.il.hr +
    sum(state == "inact" & race == "hisp") * q.inact +
    sum(state == "NODU" & race == "hisp") * q.NODU +
    sum(state == "relap" & OU.state == "preb" & race == "hisp") * q.relap.v["preb"] +
    sum(state == "relap" & OU.state == "il.lr" & race == "hisp") * q.relap.v["il.lr"] +
    sum(state == "relap" & OU.state == "il.hr" & race == "hisp") * q.relap.v["il.hr"]
  return(list(Q.TQ = Q.TQ, Q.white = Q.white, Q.black = Q.black, Q.hisp = Q.hisp)) # return the costs
}

