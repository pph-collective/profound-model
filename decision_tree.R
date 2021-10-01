#' Decision tree to determine agent state change
#' 
#' @description 
#' The decision tree stochastically determines what states the population's
#' individuals change to at a timestep. It finds the outcomes of an overdose
#' event.
#' 
#' @param od_ppl The population that has overdosed this time step.
#' @param n_nlx The number of naloxone kits available.
#' @param ou_ppl_resid Population that has overdosed by location.
#' @param params Model parameters.
#' 
#' @returns
#' outcomes determined by the model

source("naloxone_availability.R")

decision_tree <- function(od_ppl, n_nlx, ou_ppl_resid, params, seed, data) {
  list2env(params, environment())
  set.seed(seed)
  n.od <- nrow(od_ppl)
  residence <- od_ppl$residence
  col_out <- c(
    "ind", "od.death", "EMS", "hospcare", "inact", "locpriv", 
    "nlx.used", "wtns"
  )

  outcome <- matrix(0, nrow = n.od, ncol = length(col_out))
  colnames(outcome) <- col_out
  # REVIEWED ind = index; id
  outcome[, "ind"] <- od_ppl$ind
  p_nlx_avail.mx <- nlx.avail.algm(
    n_nlx,
    ou_ppl_resid,
    data$od_loc,
    data$low2priv,
    data$nlx_adj,
    data$cap
  )

  for (d in 1:n.od) {
    if (n.od == 0) {
      break
    }
    loc <- sample(c("priv", "pub"), 1, prob = data$od_loc[, residence[d]])
    locpriv <- ifelse(loc == "priv", 1, 0)
    p_wtns <- ifelse(loc == "priv", data$od_wit_priv, data$od_wit_pub)
    p_911 <- ifelse(loc == "priv", data$od_911_priv, data$od_911_pub)
    p_hosp <- data$od_hosp
    p_od2inact <- data$od_cess
    wtns <- sample.dic(p_wtns)
    p_nlx_avail <- p_nlx_avail.mx[residence[d], loc]

    if (wtns == 1) { # if witnessed
      p_nlx_avail <- sample.dic(p_nlx_avail)
      if (sample.dic(p_nlx_avail) == 1) { # naloxone, witness
        nlx.used <- 1
        EMS <- sample.dic(p_911)
        if (EMS == 1) { # EMS, naloxone, witness
          hospcare <- sample.dic(p_hosp)
          if (hospcare == 1) { # hospitalized, EMS, naloxone, witness
            od.death <- sample.dic(data$mortality_nx)
          } else { # not hospitalized, but EMS, naloxone, witness
            od.death <- sample.dic(data$mortality_nx)
          }
        } else { # not EMS, but naloxone, witness
          od.death <- sample.dic(data$mortality_nx)
          hospcare <- 0
        }
      } else { # no naloxone, but witness
        nlx.used <- 0
        EMS <- sample.dic(p_911)
        if (EMS == 1) { # EMS but no nlx
          hospcare <- sample.dic(p_hosp)
          if (hospcare == 1) { # hospitalized no nlx
            od.death <- sample.dic(data$mor_bl * data$rr_mor_EMS)
          } else { # no hospital, no nlx
            od.death <- sample.dic(data$mor_bl * data$rr_mor_EMS)
          }
        } else { # if EMS not reached (witnessed, naloxone not used by witness)
          od.death <- sample.dic(data$mor_bl)
          hospcare <- 0
        }
      }
    } else { # if not witnessed
      od.death <- sample.dic(data$mor_bl)
      EMS <- 0
      hospcare <- 0
      nlx.used <- 0
    }

    if (od.death != 1) {
      inact <- sample.dic(p_od2inact)
    } else {
      inact <- 0
    }

    outcome[d, -1] <- c(od.death, EMS, hospcare, inact, locpriv, nlx.used, wtns)
  } # end for loop

  return(outcome)
}

sample.dic <- function(prob) {
  return(sample(0:1, 1, prob = c((1 - prob), prob)))
}
