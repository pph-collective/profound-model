###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2020 #########################
###############################################################################################
# Main module for the decision tree of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: May 25, 2020
# Last update: June 01, 2020
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
# sur_bl              # baseline survival given no witness, no naloxone used, no hospital care
# rr_sur_Nx           # relative risk of survival if naloxone is used (either by wintess or EMS)
# rr_sur_EMS          # relative risk of survival with EMS if no naloxone is used
# p.od2inact          # probability to inactive (cessasion of opioid use) after surviving an overdose event

#############################################################################
# 2. Decision tree function
#############################################################################

decision.tree  <- function(od.pop, n.nlx, ou.pop.resid, seed){
  set.seed(seed)
  n.od                   <- nrow(od.pop)
  residence              <- od.pop$residence
  out.colnames           <- c("ind", "od.death", "EMS", "hospcare", "inact")
  decntree.out           <- matrix(0, nrow = n.od, ncol = length(out.colnames))
  colnames(decntree.out) <- out.colnames
  decntree.out[ , "ind"] <- od.pop$ind
  p.nlx.avail.mx         <- nlx.avail(n.nlx, ou.pop.resid, OD_loc, Low2Priv, nlx.adj)

  for (d in 1:n.od){
    loc    <- sample(c("priv", "pub"), size = 1, prob = OD_loc[ , residence[d]])
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
            od.death <- 1-sample.dic(sur_bl*rr_sur_Nx)
          } else {  # if not hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            od.death <- 1-sample.dic(sur_bl*rr_sur_Nx)
          }
        } else {  # if EMS not reached (witnessed, available, naloxone used by witness )
          od.death <- 1-sample.dic(sur_bl*rr_sur_Nx)
          hospcare <- 0
        }
      } else {  # if naloxone not used (unavailable) by witness (witnessed)
        EMS <- sample.dic(p.911)
        if (EMS == 1) {   # if EMS reached (witnessed, naloxone not used by witness)
          hospcare <- sample.dic(p.hosp)
          if (hospcare == 1) {  # if hospitalized (witnessed, naloxone not used by witness, EMS reached)
            od.death <- 1-sample.dic(sur_bl*rr_sur_EMS)
          } else {  # if not hospitalized (witnessed, naloxone not used by witness, EMS reached)
            od.death <- 1-sample.dic(sur_bl*rr_sur_EMS)
          }
        } else {  # if EMS not reached (witnessed, naloxone not used by witness)
          od.death <- 1-sample.dic(sur_bl)
          hospcare <- 0
        }
      }
    } else {  # if not witnessed
      od.death <- 1-sample.dic(sur_bl)
      EMS      <- 0
      hospcare <- 0
    }  # end if loop for decision tree
    
    if (od.death != 1){
      inact <- sample.dic(p.od2inact)
    } else {
      inact <- 0
    }
    
    decntree.out[d , -1] <- c(od.death, EMS, hospcare, inact)
  }   # end for loop
  
  return(decntree.out)
}

sample.dic <- function(prob){
  return(sample(0:1, 1, prob = c((1-prob), prob)))
}
