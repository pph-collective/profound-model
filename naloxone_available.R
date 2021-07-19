########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Function module for the microsimulation of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
#
########################################################################################
#################        Naloxone availability algorithm       #########################
########################################################################################

########################################################################################
# ref: Michael Irvine, saturation model
# The availability of naloxone in a witnessed overdose in a given setting is based on a non-linear (exponential) relationship between the number of naloxone received and the number of individuals at risk for overdose in a region 
# ou.pop.resid: number of people at risk for od in each region
# n.nlx: number of naloxone kits received by residents of each region, stratified by risk level of program
# Low2Priv: only a porportion of naloxone kits from low-risk programs go to private setting (others go to public), kits from high-risk programs all go to private (an assumption)
# cap: the maximum level naloxone availability can be in a witnssed overdose (even when reaching saturation)
# OD_loc: proportion of OD events occuring in private/public setting
# nlx.adj: adjustment term for naloxone availability, to account for occations where naloxone is in possession but not immediately accessible (or familsy/friends may fail to find)
nlx.avail.algm  <- function(n.nlx, ou.pop.resid, OD_loc, Low2Priv, nlx.adj, cap){
  n.nlx.loc             <- cbind(n.nlx["high", ] + n.nlx["low", ] * Low2Priv, n.nlx["low", ] * (1-Low2Priv))
  colnames(n.nlx.loc)   <- c("priv", "pub")
  # p.nlx.avail           <- n.nlx.loc / (ou.pop.resid$n * t(OD_loc)) * nlx.adj
  p.nlx.avail <- cap * (1 - exp(-(nlx.adj / cap) * n.nlx.loc / (ou.pop.resid$n * t(OD_loc))))
  return(p.nlx.avail)
}
