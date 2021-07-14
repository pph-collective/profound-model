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
# prob_911               # probability of seeking help from 911
# p.hosp              # probability of transporting to hospital care
# mor_bl              # baseline mortality rate (from overdose) given no witness, no naloxone used, no hospital care
# mortality_nx              # mortality if naloxone is administered
# rr_mor_EMS          # relative risk of mortality with EMS if no naloxone is used
# p.od2inact          # probability to inactive (cessasion of opioid use) after surviving an overdose event

#############################################################################
# 2. Decision tree function
#############################################################################

decisiontimesteptree <- function(od_ppl, nlx_avail, ou.ppl.resid, params, seed) {
  list2env(params, environment())
  set.seed(seed)
  num_ods <- nrow(od_ppl)
  residence <- od_ppl$residence
  out_colnames <- c("id", "death", "EMS", "hospcare", "inact", "private_loc", "nlx_used", "wtns")
  decntree.out <- matrix(0, nrow = num_ods, ncol = length(out_colnames))
  colnames(decntree.out) <- out_colnames

  decntree.out[, "id"] <- od_ppl$id
  p.nlx.avail.mx <- nlx.avail.algm(nlx_avail, ou.ppl.resid, OD_loc, Low2Priv, nlx.adj, cap)

  for (d in 1:num_ods) {
    loc <- sample(c("priv", "pub"), size = 1, prob = OD_loc[, residence[d]])
    private_loc <- ifelse(loc == "priv", 1, 0)
    p.wtns <- ifelse(loc == "priv", OD_wit_priv, OD_wit_pub)
    prob_911 <- ifelse(loc == "priv", OD_911_priv, OD_911_pub)
    p.hosp <- OD_hosp
    p.od2inact <- OD_cess

    p.nlx.avail <- p.nlx.avail.mx[residence[d], loc]

    wtns <- occurs(p.wtns)
    if (wtns == 1) { # if witnessed
      if (occurs(p.nlx.avail) == 1) { # if naloxone used by witness
        nlx_used <- 1
        EMS <- occurs(prob_911)
        if (EMS == 1) { # if EMS reached (witnessed, available, naloxone used by witness)
          hospcare <- occurs(p.hosp)
          if (hospcare == 1) { # if hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            death <- occurs(mortality_nx)
          } else { # if not hospitalized (witnessed, available, naloxone used by witness, EMS reached)
            death <- occurs(mortality_nx)
          }
        } else { # if EMS not reached (witnessed, available, naloxone used by witness)
          death <- occurs(mortality_nx)
          hospcare <- 0
        }
      } else { # if naloxone not used (unavailable) by witness (witnessed)
        nlx_used <- 0
        EMS <- occurs(prob_911)
        if (EMS == 1) { # if EMS reached (witnessed, naloxone not used by witness)
          hospcare <- occurs(p.hosp)
          if (hospcare == 1) { # if hospitalized (witnessed, naloxone not used by witness, EMS reached)
            death <- occurs(mor_bl * rr_mor_EMS)
          } else { # if not hospitalized (witnessed, naloxone not used by witness, EMS reached)
            death <- occurs(mor_bl * rr_mor_EMS)
          }
        } else { # if EMS not reached (witnessed, naloxone not used by witness)
          death <- occurs(mor_bl)
          hospcare <- 0
        }
      }
    } else { # if not witnessed
      death <- occurs(mor_bl)
      EMS <- 0
      hospcare <- 0
      nlx_used <- 0
    } # end if loop for decision tree

    if (death != 1) {
      inact <- occurs(p.od2inact)
    } else {
      inact <- 0
    }

    decntree.out[d, -1] <- c(death, EMS, hospcare, inact, private_loc, nlx_used, wtns)
  } # end for loop

  return(decntree.out)
}

occurs <- function(prob) {
  return(sample(0:1, 1, prob = c((1 - prob), prob)))
}
