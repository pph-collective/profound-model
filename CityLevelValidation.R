ODdeaths16 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
ODdeaths17 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
ODdeaths18 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
ODdeaths19 <- matrix(0, nrow = n.rgn, ncol = nrow(calib.result.mx))
row.names(ODdeaths16) <- row.names(ODdeaths17) <- row.names(ODdeaths18) <- row.names(ODdeaths19) <- v.rgn
for (ss in 1:nrow(calib.result.mx)){
  print(paste0("Parameter set: ", ss))
  calib.seed = calib.result.mx[ss, "seed"]

  for (pp in 1:length(nm.calp)){
    assign(nm.calp[pp], calib.result.mx[ss,nm.calp[pp]])
  }
  ## Fentanyl use status for initial population determined externally (allow to vary) ##
  n.opioid <- sum(init.pop$curr.state != "NODU")
  n.noud   <- sum(init.pop$curr.state == "NODU")
  # determine fentanyl use among initial population who use opioids 
  set.seed(calib.seed)
  fx         <- sample(0:1, size = n.opioid, prob = c(1-ini.OUD.fx, ini.OUD.fx), replace = T)
  init.pop$fx[init.pop$curr.state != "NODU"] <- fx
  # determine fentanyl use among initial population who use stimulants (non-opioid)
  set.seed(calib.seed*2)
  fx         <- sample(0:1, size = n.noud, prob = c(1-ini.NOUD.fx, ini.NOUD.fx), replace = T)
  init.pop$fx[init.pop$curr.state == "NODU"] <- fx
  
  # Overdose probability matrix (per month)
  od.matrix             <- matrix(0, nrow = 4, ncol = 2)
  rownames(od.matrix)   <- c("preb", "il.lr", "il.hr", "NODU")
  colnames(od.matrix)   <- c("first", "subs")
  od.matrix["preb", "subs"]   <- od.preb.sub
  od.matrix["il.lr", "subs"]  <- od.il.lr.sub
  od.matrix["il.hr", "subs"]  <- od.il.lr.sub * multi.hr
  od.matrix["NODU", "subs"]   <- od.NODU.sub
  od.matrix[ , "first"]       <- od.matrix[ , "subs"] / multi.sub
  
  sim_sq    <- MicroSim(init.pop, n.t, v.state, d.c, PT.out = TRUE, Str = "SQ", seed = calib.seed)        # run for no treatment

  ODdeaths16[ , ss] <- colSums(sim_sq$m.oddeath[1:12, ])
  ODdeaths17[ , ss] <- colSums(sim_sq$m.oddeath[13:24, ])
  ODdeaths18[ , ss] <- colSums(sim_sq$m.oddeath[25:36, ])
  ODdeaths19[ , ss] <- colSums(sim_sq$m.oddeath[37:48, ])
}

write.xlsx(ODdeaths16, file="CityLevelValidatin.xlsx", sheetName = "2016",
           col.names = F, row.names = T)
write.xlsx(ODdeaths17, file="CityLevelValidatin.xlsx", sheetName = "2017", append=TRUE,
           col.names = F, row.names = T)
write.xlsx(ODdeaths18, file="CityLevelValidatin.xlsx", sheetName = "2018", append=TRUE,
           col.names = F, row.names = T)
write.xlsx(ODdeaths19, file="CityLevelValidatin.xlsx", sheetName = "2019", append=TRUE,
           col.names = F, row.names = T)