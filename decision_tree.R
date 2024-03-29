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
#########################           Decision Tree           ####################################
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
# REVIEWED if no witness, no naloxone \ can we just say probability of naloxone in the variable name (prob_nlx)?
# p.nlx.wtns          # probability naloxone is used/administered by witness if available
# p.911               # probability of seeking help from 911
# p.hosp              # probability of transporting to hospital care
# mor_bl              # baseline mortality rate (from overdose) given no witness, no naloxone used, no hospital care
# mor_nx              # mortality if naloxone is administered
# rr_mor_EMS          # relative risk of mortality with EMS if no naloxone is used
# p.od2inact          # probability to inactive (cessasion of opioid use) after surviving an overdose event

#############################################################################
# 2. Decision tree function
#############################################################################

decision_tree <- function(od_ppl, p.nlx.avail.mx, params, strategy = "SQ", t, seed) {
  # REVIEWED n.nlx is naloxone available in city; resid is study ppl in city, change params to params or parameters
  list2env(params, environment())
  set.seed(seed)
  n.od <- nrow(od_ppl)
  residence <- od_ppl$residence
  OU.state <- od_ppl$OU.state
  out.colnames <- c("ind", "od.death", "EMS", "hospcare", "inact", "locpriv", "nlx.used", "wtns")
  decntree.out <- matrix(0, nrow = n.od, ncol = length(out.colnames))
  colnames(decntree.out) <- c("ind", "od.death", "EMS", "hospcare", "inact", "locpriv", "nlx.used", "wtns")
  # REVIEWED ind = index; id
  decntree.out[, "ind"] <- od_ppl$ind

  for (d in 1:n.od) {
    loc <- sample(c("priv", "pub"), size = 1, prob = c(1-OD_loc_pub, OD_loc_pub))
    locpriv <- ifelse(loc == "priv", 1, 0)
    p.wtns <- ifelse(loc == "priv", OD_wit_priv, OD_wit_pub)
    p.911 <- ifelse(loc == "priv", OD_911_priv, OD_911_pub)
    p.hosp <- OD_hosp
    p.od2inact <- OD_cess
    if (!grepl("10K", strategy) | t<=48){
      p.nlx.avail <- as.numeric(p.nlx.avail.mx[residence[d]])
    } else {
      p.nlx.avail <- p.nlx.avail.mx[residence[d], OU.state[d]]
    }
    wtns <- sample.dic(p.wtns)
    if (wtns == 1) { # if witnessed
      nlx.used <- sample.dic(p.nlx.avail)
      if (nlx.used == 1) { # if naloxone available by witness (witnessed)
        EMS <- sample.dic(p.911)
        if (EMS == 1) { # if EMS reached (witnessed, available, naloxone used by witness )
          hospcare <- sample.dic(p.hosp)
          if (hospcare == 1) { # if hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            od.death <- sample.dic(mor_nx)
          } else { # if not hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            od.death <- sample.dic(mor_nx)
          }
        } else { # if EMS not reached (witnessed, available, naloxone used by witness )
          od.death <- sample.dic(mor_nx)
          hospcare <- 0
        }
      } else { # if naloxone not used (unavailable) by witness (witnessed)
        nlx.used <- 0
        EMS <- sample.dic(p.911)
        if (EMS == 1) { # if EMS reached (witnessed, naloxone not used by witness)
          hospcare <- sample.dic(p.hosp)
          if (hospcare == 1) { # if hospitalized (witnessed, naloxone not used by witness, EMS reached)
            od.death <- sample.dic(mor_bl * rr_mor_EMS)
          } else { # if not hospitalized (witnessed, naloxone not used by witness, EMS reached)
            od.death <- sample.dic(mor_bl * rr_mor_EMS)
          }
        } else { # if EMS not reached (witnessed, naloxone not used by witness)
          od.death <- sample.dic(mor_bl)
          hospcare <- 0
        }
      }
    } else { # if not witnessed
      od.death <- sample.dic(mor_bl)
      EMS <- 0
      hospcare <- 0
      nlx.used <- 0
    } # end if loop for decision tree

    if (od.death != 1 & od_ppl$OU.state[d] != "NODU" & hospcare == 1) {
      inact <- sample.dic(p.od2inact)
    } else {
      inact <- 0
    }

    decntree.out[d, -1] <- c(od.death, EMS, hospcare, inact, locpriv, nlx.used, wtns)
  } # end for loop

  return(decntree.out)
}

sample.dic <- function(prob) {
  return(sample(0:1, 1, prob = c((1 - prob), prob)))
}
