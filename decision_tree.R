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
  n_od <- nrow(od_ppl)
  residence <- od_ppl$residence
  col_out <- c(
    "ind", "od_death", "EMS", "hospcare", "inact", "locpriv",
    "nlx_used", "wtns"
  )

  outcome <- matrix(0, nrow = n_od, ncol = length(col_out))
  colnames(outcome) <- col_out
  outcome[, "ind"] <- od_ppl$ind
  p_nlx_avail_mx <- nlx_avail_algm(
    n_nlx,
    ou_ppl_resid,
    data$od_loc,
    data$low2priv,
    data$nlx_adj,
    data$cap
  )

  for (d in 1:n_od) {
    if (n_od == 0) {
      break
    }
    loc <- sample(c("priv", "pub"), 1, prob = data$od_loc[, residence[d]])
    locpriv <- ifelse(loc == "priv", 1, 0)
    p_wtns <- ifelse(loc == "priv", data$od_wit_priv, data$od_wit_pub)
    p_911 <- ifelse(loc == "priv", data$od_911_priv, data$od_911_pub)
    p_hosp <- data$od_hosp
    p_od2inact <- data$od_cess
    wtns <- sample_dic(p_wtns)
    p_nlx_avail <- p_nlx_avail_mx[residence[d], loc]

    if (wtns == 1) { # if witnessed
      p_nlx_avail <- sample_dic(p_nlx_avail)
      if (sample_dic(p_nlx_avail) == 1) { # naloxone, witness
        nlx_used <- 1
        ems <- sample_dic(p_911)
        if (ems == 1) { # EMS, naloxone, witness
          hospcare <- sample_dic(p_hosp)
          if (hospcare == 1) { # hospitalized, EMS, naloxone, witness
            od_death <- sample_dic(data$mortality_nx)
          } else { # not hospitalized, but EMS, naloxone, witness
            od_death <- sample_dic(data$mortality_nx)
          }
        } else { # not EMS, but naloxone, witness
          od_death <- sample_dic(data$mortality_nx)
          hospcare <- 0
        }
      } else { # no naloxone, but witness
        nlx_used <- 0
        ems <- sample_dic(p_911)
        if (ems == 1) { # EMS but no nlx
          hospcare <- sample_dic(p_hosp)
          if (hospcare == 1) { # hospitalized no nlx
            od_death <- sample_dic(data$mor_bl * data$rr_mor_ems)
          } else { # no hospital, no nlx
            od_death <- sample_dic(data$mor_bl * data$rr_mor_ems)
          }
        } else { # if EMS not reached (witnessed, naloxone not used by witness)
          od_death <- sample_dic(data$mor_bl)
          hospcare <- 0
        }
      }
    } else { # if not witnessed
      od_death <- sample_dic(data$mor_bl)
      ems <- 0
      hospcare <- 0
      nlx_used <- 0
    }

    if (od_death != 1) {
      inact <- sample_dic(p_od2inact)
    } else {
      inact <- 0
    }

    outcome[d, -1] <- c(od_death, ems, hospcare, inact, locpriv, nlx_used, wtns)
  } # end for loop
  return(outcome)
}

sample_dic <- function(prob) {
  return(sample(0:1, 1, prob = c((1 - prob), prob)))
}
