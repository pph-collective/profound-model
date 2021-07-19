# TODO file header

ODdeaths16 <- matrix(0, nrow = num_regions, ncol = nrow(calib.result.mx))
ODdeaths17 <- matrix(0, nrow = num_regions, ncol = nrow(calib.result.mx))
ODdeaths18 <- matrix(0, nrow = num_regions, ncol = nrow(calib.result.mx))
ODdeaths19 <- matrix(0, nrow = num_regions, ncol = nrow(calib.result.mx))
row.names(ODdeaths16) <- row.names(ODdeaths17) <- row.names(ODdeaths18) <- row.names(ODdeaths19) <- v.region
for (ss in 1:nrow(calib.result.mx)) {
  print(paste0("Parameter set: ", ss))
  calib.seed <- calib.result.mx[ss, "seed"]

  for (pp in 1:length(cal_param_names)) {
    assign(cal_param_names[pp], calib.result.mx[ss, cal_param_names[pp]])
  }
  ## Fentanyl use status for initial population determined externally (allow to vary) ##
  # TO_REVIEW is there ever a population variable like this that isn't the initial population? Can we change the variable name to just pop (ppl) instead of init_ppl (init_ppl)
  num_opioid <- sum(init_ppl$curr.state != "NODU")
  n.noud <- sum(init_ppl$curr.state == "NODU")
  # determine fentanyl use among initial population who use opioids
  set.seed(calib.seed)
  fx <- sample(0:1, size = num_opioid, prob = c(1 - ini.OUD.fx, ini.OUD.fx), replace = T)
  init_ppl$fx[init_ppl$curr.state != "NODU"] <- fx
  # determine fentanyl use among initial population who use stimulants (non-opioid)
  set.seed(calib.seed * 2)
  fx <- sample(0:1, size = n.noud, prob = c(1 - ini.NOUD.fx, ini.NOUD.fx), replace = T)
  init_ppl$fx[init_ppl$curr.state == "NODU"] <- fx

  # Overdose probability matrix (per month)
  # TO_REVIEW sub/subs?
  overdose_probs <- matrix(0, nrow = 4, ncol = 2)
  rownames(overdose_probs) <- c("preb", "il.lr", "il.hr", "NODU")
  colnames(overdose_probs) <- c("first", "subs")
  overdose_probs["preb", "subs"] <- od.preb.sub
  overdose_probs["il.lr", "subs"] <- od.il.lr.sub
  overdose_probs["il.hr", "subs"] <- od.il.lr.sub * multi.hr
  overdose_probs["NODU", "subs"] <- od.NODU.sub
  overdose_probs[, "first"] <- overdose_probs[, "subs"] / multi.sub
  # TO_REVIEW sim_sq status quo?
  # run status quo simulation
  sim_sq <- MicroSim(init_ppl, timesteps, agent_states, d.c, PT.out = TRUE, strategy = "SQ", seed = calib.seed) # run for no treatment

  ODdeaths16[, ss] <- colSums(sim_sq$m.oddeath[1:12, ])
  ODdeaths17[, ss] <- colSums(sim_sq$m.oddeath[13:24, ])
  ODdeaths18[, ss] <- colSums(sim_sq$m.oddeath[25:36, ])
  ODdeaths19[, ss] <- colSums(sim_sq$m.oddeath[37:48, ])
}

write.xlsx(ODdeaths16,
  file = "CityLevelValidation.xlsx", sheetName = "2016",
  col.names = F, row.names = T
)
write.xlsx(ODdeaths17,
  file = "CityLevelValidation.xlsx", sheetName = "2017", append = TRUE,
  col.names = F, row.names = T
)
write.xlsx(ODdeaths18,
  file = "CityLevelValidation.xlsx", sheetName = "2018", append = TRUE,
  col.names = F, row.names = T
)
write.xlsx(ODdeaths19,
  file = "CityLevelValidation.xlsx", sheetName = "2019", append = TRUE,
  col.names = F, row.names = T
)
