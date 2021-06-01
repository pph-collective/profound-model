###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2020 #########################
###############################################################################################
# Main module for the decision tree of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: May 25, 2020
# Last update: April 12, 2021
#
###############################################################################################
#########################           Decsion Tree           ####################################
###############################################################################################

###############################################################################################
####    Decision tree to determine consequence following overdose                          ####
####    6 health states: prescribed, illicit (L/H), inactive, non-opioid, relapsed, death  ####
####    1 health event:  Overdose                                                          ####
####    Attributes:      state, age, sex, fentanyl, overdosed, pre.state,                  ####
####    Built to inform Naloxone distribution strategies to prevent overdsoe death         ####
###############################################################################################

#############################################################################
# 1. Decision tree parameters
#############################################################################
# INPUT PARAMETERS
# p.wtns              # probability an overodse is witnessed
# p.nlx.wtns          # probability naloxone is used/administered by witness if available             
# p.911               # probability of seeking help from 911
# p.hosp              # probability of transporting to hospital care
# mor_bl              # baseline mortality rate (from overdose) given no witness, no naloxone used, no hospital care
# mor_Nx              # mortality if naloxone is administered
# rr_mor_EMS          # relative risk of mortality with EMS if no naloxone is used
# p.od2inact          # probability to inactive (cessasion of opioid use) after surviving an overdose event

#############################################################################
# 2. Decision tree function
#############################################################################

decision.tree  <- function(od.pop, n.nlx, ou.pop.resid, vparameters, seed){
  list2env(vparameters, environment())
  set.seed(seed)
  n.od                   <- nrow(od.pop)
  residence              <- od.pop$residence
  out.colnames           <- c("ind", "od.death", "EMS", "hospcare", "inact", "locpriv")
  decntree.out           <- matrix(0, nrow = n.od, ncol = length(out.colnames))
  colnames(decntree.out) <- out.colnames
  decntree.out[ , "ind"] <- od.pop$ind
  p.nlx.avail.mx         <- nlx.avail(n.nlx, ou.pop.resid, OD_loc, Low2Priv, nlx.adj, cap)

  for (d in 1:n.od){
    loc    <- sample(c("priv", "pub"), size = 1, prob = OD_loc[ , residence[d]])
    locpriv<- ifelse(loc == "priv", 1, 0)
    p.wtns <- ifelse(loc == "priv", OD_wit_priv, OD_wit_pub)
    p.911  <- ifelse(loc == "priv", OD_911_priv, OD_911_pub)
    p.hosp <- OD_hosp
    p.od2inact <- OD_cess
    
    wtns   <- sample.dic(p.wtns)
    p.nlx.avail <- p.nlx.avail.mx[residence[d], loc]
    
    if (wtns == 1) {   # if witnessed
      nlx.avail <- sample.dic(p.nlx.avail)
      if (nlx.avail == 1){   # if naloxone available by witness (witnessed)
        EMS <- sample.dic(p.911)
        if (EMS == 1) {  # if EMS reached (witnessed, available, naloxone used by witness )
          hospcare <- sample.dic(p.hosp)
          if (hospcare == 1) {  # if hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            od.death <- sample.dic(mor_Nx)
          } else {  # if not hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            od.death <- sample.dic(mor_Nx)
          }
        } else {  # if EMS not reached (witnessed, available, naloxone used by witness )
          od.death <- sample.dic(mor_Nx)
          hospcare <- 0
        }
      } else {  # if naloxone not used (unavailable) by witness (witnessed)
        EMS <- sample.dic(p.911)
        if (EMS == 1) {   # if EMS reached (witnessed, naloxone not used by witness)
          hospcare <- sample.dic(p.hosp)
          if (hospcare == 1) {  # if hospitalized (witnessed, naloxone not used by witness, EMS reached)
            od.death <- sample.dic(mor_bl*rr_mor_EMS)
          } else {  # if not hospitalized (witnessed, naloxone not used by witness, EMS reached)
            od.death <- sample.dic(mor_bl*rr_mor_EMS)
          }
        } else {  # if EMS not reached (witnessed, naloxone not used by witness)
          od.death <- sample.dic(mor_bl)
          hospcare <- 0
        }
      }
    } else {  # if not witnessed
      od.death <- sample.dic(mor_bl)
      EMS      <- 0
      hospcare <- 0
    }  # end if loop for decision tree
    
    if (od.death != 1){
      inact <- sample.dic(p.od2inact)
    } else {
      inact <- 0
    }
    
    decntree.out[d , -1] <- c(od.death, EMS, hospcare, inact, locpriv)
  }   # end for loop
  
  return(decntree.out)
}

sample.dic <- function(prob){
  return(sample(0:1, 1, prob = c((1-prob), prob)))
}
