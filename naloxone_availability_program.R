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
nlx.avail.algm <- function(n.nlx, crosstab_drug_resid, OD_loc_pub, nlx.adj, cap, eff.pharNlx, strategy = "SQ", t) {
  if (!grepl("10K", strategy)){
    n.nlx.total <- n.nlx$OEND + n.nlx$Pharm / eff.pharNlx
    p.nlx.avail <- cap * (1 - exp(-(nlx.adj / cap) * n.nlx.total / rowSums(crosstab_drug_resid)))
    return(p.nlx.avail)
  } else if(t <= 48){
    n.nlx$OEND <- n.nlx$OEND[-length(n.nlx$OEND)]
    n.nlx.total <- n.nlx$OEND + n.nlx$Pharm / eff.pharNlx
    p.nlx.avail <- cap * (1 - exp(-(nlx.adj / cap) * n.nlx.total / rowSums(crosstab_drug_resid)))
    return(p.nlx.avail)
  } else {
    tenk <- sum(n.nlx$OEND[length(n.nlx$OEND)])
    n.nlx.base <- n.nlx$OEND[-length(n.nlx$OEND)] + n.nlx$Pharm / eff.pharNlx
    n.nlx.drug.resid.base <- as.vector(t(n.nlx.base)) * crosstab_drug_resid / rowSums(crosstab_drug_resid)
    if (strategy == "SSP_10K"){
      SSP_drug_resid <- crosstab_drug_resid[, c("il.hr", "NODU")]*rep(c(1, 0.133), each = nrow(crosstab_drug_resid))
      n.nlx.drug.resid <- n.nlx.drug.resid.base
      n.nlx.drug.resid[, c("il.hr", "NODU")] <- tenk * programprop$SSP * SSP_drug_resid / rowSums(SSP_drug_resid) + n.nlx.drug.resid.base[, c("il.hr", "NODU")]
    } else if (strategy == "Outreach_10K"){
      Outreach_drug_resid <- crosstab_drug_resid[, c("il.lr", "il.hr", "NODU")]
      n.nlx.drug.resid <- n.nlx.drug.resid.base
      n.nlx.drug.resid[, c("il.lr", "il.hr", "NODU")] <- tenk * programprop$Outreach * Outreach_drug_resid / rowSums(Outreach_drug_resid) + n.nlx.drug.resid.base[, c("il.lr", "il.hr", "NODU")]
    } else if (strategy == "MailEvent_10K"){
      n.nlx.drug.resid <- n.nlx.drug.resid.base
      n.nlx.drug.resid <- tenk * programprop$MailEvent * crosstab_drug_resid / rowSums(crosstab_drug_resid) + n.nlx.drug.resid.base
    } else if (strategy == "Healthcare_10K"){
      Healthcare_drug_resid <- crosstab_drug_resid[, c("preb")]
      n.nlx.drug.resid <- n.nlx.drug.resid.base
      n.nlx.drug.resid[, c("preb")] <- tenk * programprop$Healthcare  + n.nlx.drug.resid.base[, c("preb")]
    }
    p.nlx.avail <- cap * (1 - exp(-(nlx.adj / cap) * n.nlx.drug.resid / crosstab_drug_resid))
    return(p.nlx.avail)
  }
}
