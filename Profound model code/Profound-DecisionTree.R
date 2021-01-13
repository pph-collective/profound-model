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
# 1. SET directpry and workspace
#############################################################################

# INPUT PARAMETERS
p.wtns              <- 0.75           # probability an overodse is witnessed
p.nlx.avail.bl      <- 0.4            # probability naloxone is available during a witnessed overdose, time-varying and stratified by region, will be updated in decision tree
p.nlx.wtns          <- 0.9            # probability naloxone is used/administered by witness if available             
p.ems.nlx           <- 0.85           # probability of seeking help from EMS when naloxone is used
p.ems.Nnlx          <- 0.6            # probability of seeking help from EMS when naloxone is not used (N stands for "no")
p.nlx.ems           <- 0.89           # probability naloxone is used/administered by EMS
p.hospcare          <- 0.9            # probability of transporting to hospital care
mor.N               <- 1-0.88         # morality given no witness, no naloxone used, no hospital care (N stands for "no")
rr.mor.nlx          <- 0.48           # relative risk of mortality if naloxone is used (either by wintess or EMS)
rr.mor.ems          <- 0.8            # relative risk of mortality with EMS if no naloxone is used
p.od2inact          <- 0.15           # probability to inactive (cessasion of opioid use) after surviving an overdose event
rr.inact.hospcare   <- 1.2            # relative risk to inactive (cessasion of opioid use) if received hospital care

#############################################################################
# 2. Decision tree function
#############################################################################

decision.tree  <- function(od.pop, n.nlx.v, n.od_death.v, acc.nlx.matrix, seed){
  set.seed(seed)
  n.od                   <- nrow(od.pop)
  residence              <- od.pop$residence
  out.colnames           <- c("ind", "od.death", "EMS", "hospcare", "inact")
  decntree.out           <- matrix(0, nrow = n.od, ncol = length(out.colnames))
  colnames(decntree.out) <- out.colnames
  decntree.out[ , "ind"] <- od.pop$ind
  p.nlx.avail.v          <- nlx.avail(n.nlx.v, n.od_death.v, acc.nlx.matrix)
  
  for (d in 1:n.od){
    wtns <- sample.dic(p.wtns)
    if (wtns == 1) {   # if witnessed
      p.nlx.avail <- p.nlx.avail.v[residence[d] == v.rgn]
      nlx.avail <- sample.dic(p.nlx.avail)
      if (nlx.avail == 1){   # if naloxone available by witness (witnessed)
        nlx.wtns <- sample.dic(p.nlx.wtns)
        if (nlx.wtns == 1) {  # if naloxone used by witness (witnessed, available)
          EMS <- sample.dic(p.ems.nlx)
          if (EMS == 1) {  # if EMS reached (witnessed, available, naloxone used by witness )
            hospcare <- sample.dic(p.hospcare)
            if (hospcare == 1) {  # if hospitalized (witnessed, available, naloxone used by witness, EMS reached)
              od.death <- sample.dic(mor.N*rr.mor.nlx)
            } else {  # if not hospitalized (witnessed, available, naloxone used by witness, EMS reached)
              od.death <- sample.dic(mor.N*rr.mor.nlx)
            }
          } else {  # if EMS not reached (witnessed, available, naloxone used by witness )
            od.death <- sample.dic(mor.N*rr.mor.nlx)
            hospcare <- 0
          }
        } else {  # if naloxone not used by witness (witnessed, available)
          EMS <- sample.dic(p.ems.Nnlx)
          if (EMS == 1) {  # if EMS reached (witnessed, available, naloxone not used by witness)
            nlx.ems  <- sample.dic(p.nlx.ems)
            if (nlx.ems == 1) {  # if naloxone used by EMS (witnessed, available, naloxone not used by witness, EMS reached)
              hospcare <- sample.dic(p.hospcare)
              if (hospcare == 1) {  # if hospitalized (witnessed, available, naloxone not used by witness, EMS reached, naloxone used by EMS )
                od.death <- sample.dic(mor.N*rr.mor.nlx)
              } else {  # if not hospitalized (witnessed, available, naloxone not used by witness, EMS reached, naloxone used by EMS )
                od.death <- sample.dic(mor.N*rr.mor.nlx)
              }
            } else {  # if naloxone not used by EMS (witnessed, available, naloxone not used by witness, EMS reached)
              hospcare <- sample.dic(p.hospcare)
              if (hospcare == 1) {  # if hospitalized(witnessed, available, naloxone not used by witness, EMS reached, naloxone not used by EMS)
                od.death <- sample.dic(mor.N*rr.mor.ems)
              } else {  # if not hospitalized(witnessed, available, naloxone not used by witness, EMS reached, naloxone not used by EMS)
                od.death <- sample.dic(mor.N*rr.mor.ems)
              }
            }
          } else {  # if EMS not reached (witnessed, available, naloxone not used by witness)
            od.death <- sample.dic(mor.N)
            hospcare <- 0
          }
        }
      } else {  # if naloxone not used (unavailable) by witness (witnessed)
        EMS <- sample.dic(p.ems.Nnlx)
        if (EMS == 1) {   # if EMS reached (witnessed, naloxone not used by witness)
          nlx.ems  <- sample.dic(p.nlx.ems)
          if (nlx.ems == 1) {  # if naloxone used by EMS (witnssed, naloxone not used by witness, EMS reached)
            hospcare <- sample.dic(p.hospcare)
            if (hospcare == 1) {  # if hospitalized (witnessed, naloxone not used by witness, EMS reached, naloxone used by EMS)
              od.death <- sample.dic(mor.N*rr.mor.nlx)
            } else {  # if not hospitalized (witnessed, naloxone not used by witness, EMS reached, naloxone used by EMS)
              od.death <- sample.dic(mor.N*rr.mor.nlx)
            }
          } else {  # if naloxone not used by EMS (witnssed, naloxone not used by witness, EMS reached)
            hospcare <- sample.dic(p.hospcare)
            if (hospcare == 1) {  # if hospitalized (witnessed, naloxone not used by witness, EMS reached, naloxone not used by EMS)
              od.death <- sample.dic(mor.N*rr.mor.ems)
            } else {  # if not hospitalized (witnessed, naloxone not used by witness, EMS reached, naloxone not used by EMS)
              od.death <- sample.dic(mor.N*rr.mor.ems)
            }
          }
        } else {  # if EMS not reached (witnessed, naloxone not used by witness)
          od.death <- sample.dic(mor.N)
          hospcare <- 0
        }
      }
    } else {  # if not witnessed
      od.death <- sample.dic(mor.N)
      EMS      <- 0
      hospcare <- 0
    }  # end if loop for decision tree
    
    if (od.death != 1){
      if (hospcare != 1){
        inact <- sample.dic(p.od2inact)
      } else {
        inact <- sample.dic(p.od2inact*rr.inact.hospcare)
      }
    } else {
      inact <- 0
    }
    
    decntree.out[d , -1] <- c(od.death, EMS, hospcare, inact)
  }   # end for loop
  
  return(decntree.out)
}

nlx.avail      <- function(n.nlx.v, n.od_death.v, acc.nlx.matrix){
  n.nlx.to              <- n.nlx.v %*% acc.nlx.matrix
  p.nlx.avail           <- p.nlx.avail.bl * n.nlx.to / n.od_death.v * nlx.adj
  colnames(p.nlx.avail) <- v.rgn
  return(p.nlx.avail)
}

sample.dic <- function(prob){
  return(sample(0:1, 1, prob = c((1-prob), prob)))
}
