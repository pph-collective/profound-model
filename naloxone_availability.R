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
# low2priv: only a porportion of naloxone kits from low-risk programs go to private setting (others go to public), kits from high-risk programs all go to private (an assumption)
# cap: the maximum level naloxone availability can be in a witnssed overdose (even when reaching saturation)
# od_loc: proportion of OD events occuring in private/public setting
# nlx_adj: adjustment term for naloxone availability, to account for occations where naloxone is in possession but not immediately accessible (or familsy/friends may fail to find)
nlx_avail_algm <- function(n.nlx, ou.pop.resid, od_loc, low2priv, nlx_adj, cap) {
  n.nlx.loc <- cbind(n.nlx["high", ] + n.nlx["low", ] * low2priv, n.nlx["low", ] * (1 - low2priv))
  colnames(n.nlx.loc) <- c("priv", "pub")
  # ou <- merge(ou.pop.resid, t(od_loc), by.x = "residence", by.y = 0, all = TRUE)
  # p_nlx_avail <- cap * (1 - exp(-(nlx_adj / cap) * n.nlx.loc / (ou$n * t(od_loc))))
  p.nlx.avail <- cap * (1 - exp(-(nlx_adj / cap) * n.nlx.loc / (ou.pop.resid$n * t(od_loc))))
  # return(p_nlx_avail)
}
