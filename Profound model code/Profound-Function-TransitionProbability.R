########################################################################################
################# PROFOUND Naloxone Distribution model #### 2020 #######################
########################################################################################
# Function module for the microsimulation of the Profound Naloxone distribution model: 
#
# Author: Xiao Zang, PhD; Shayla Nolen, MPH
# Marshall Lab, Department of Epidemiology, Brown University
#
# Created: May 16, 2020
# Last update: May 30, 2020
#
########################################################################################
#################        Transistion probability function      #########################
########################################################################################

trans.prob <- function(pop.t){
  # initialize a matrix to store transition probabilities of each individual
  n.state <- n.state
  trans.prob.matrix <- matrix(NA, n.state + 1, nrow(pop.t)) #add one row for OD
  rownames(trans.prob.matrix) <- c(v.state, "od")
  
  # create a vector to store baseline mortality (excluding od death) for each individual according to age and treatment
  mor.rate <- numeric(nrow(pop.t))                                 
  mor.rate[filter(pop.t, age <= 35 & curr.state != "inact")$ind] <- mor.matrix["untreated", "ag1"]
  mor.rate[filter(pop.t, age <= 35 & curr.state == "inact")$ind] <- mor.matrix["treated",   "ag1"]
  mor.rate[filter(pop.t, age > 35  & curr.state != "inact")$ind] <- mor.matrix["untreated", "ag2"]
  mor.rate[filter(pop.t, age > 35  & curr.state == "inact")$ind] <- mor.matrix["treated",   "ag2"]
  
  # create a vector to store probability of overdose for each individual according to ever overdosed and fentanyl
  od.rate <- numeric(nrow(pop.t))                                 
  od.rate[filter(pop.t, curr.state == "preb"  & ever.od == 0 & fx == 0)$ind]  <- od.matrix["preb", "first"]
  od.rate[filter(pop.t, curr.state == "preb"  & ever.od == 1 & fx == 0)$ind]  <- od.matrix["preb", "subs"]
  od.rate[filter(pop.t, curr.state == "preb"  &                fx == 1)$ind]  <- od.matrix["preb", "fx"]
  od.rate[filter(pop.t, curr.state == "il.lr" & ever.od == 0 & fx == 0)$ind]  <- od.matrix["il.lr", "first"]
  od.rate[filter(pop.t, curr.state == "il.lr" & ever.od == 1 & fx == 0)$ind]  <- od.matrix["il.lr", "subs"]
  od.rate[filter(pop.t, curr.state == "il.lr" &                fx == 1)$ind]  <- od.matrix["il.lr", "fx"]
  od.rate[filter(pop.t, curr.state == "il.hr" & ever.od == 0 & fx == 0)$ind]  <- od.matrix["il.hr", "first"]
  od.rate[filter(pop.t, curr.state == "il.hr" & ever.od == 1 & fx == 0)$ind]  <- od.matrix["il.hr", "subs"]
  od.rate[filter(pop.t, curr.state == "il.hr" &                fx == 1)$ind]  <- od.matrix["il.hr", "fx"]
  od.rate[filter(pop.t, curr.state == "NODU"  & ever.od == 0 & fx == 0)$ind]  <- od.matrix["NODU", "first"]
  od.rate[filter(pop.t, curr.state == "NODU"  & ever.od == 1 & fx == 0)$ind]  <- od.matrix["NODU", "subs"]
  od.rate[filter(pop.t, curr.state == "NODU"  &                fx == 1)$ind]  <- od.matrix["NODU", "fx"]
  od.rate[filter(pop.t, curr.state == "relap" & OU.state == "preb")$ind]      <- od.matrix["preb",  "relap"]
  od.rate[filter(pop.t, curr.state == "relap" & OU.state == "il.lr")$ind]     <- od.matrix["il.lr", "relap"]
  od.rate[filter(pop.t, curr.state == "relap" & OU.state == "il.hr")$ind]     <- od.matrix["il.hr", "relap"]
  od.rate[filter(pop.t, curr.state == "relap" & OU.state == "NODU")$ind]      <- od.matrix["NODU",  "relap"]
  
  # update the trans.prob matrix with the appropriate probabilities
  ind.preb <- pop.t$curr.state == "preb"
  if (sum(ind.preb) != 0){
    trans.prob.matrix[ , ind.preb]   <- rbind(1 - p.preb2il.lr - p.preb2inact - mor.rate[ind.preb] - od.rate[ind.preb], 
                                              rep(p.preb2il.lr, sum(ind.preb)),
                                              rep(0, sum(ind.preb)),
                                              rep(p.preb2inact, sum(ind.preb)),
                                              rep(0, sum(ind.preb)), 
                                              rep(0, sum(ind.preb)), 
                                              mor.rate[ind.preb],
                                              od.rate[ind.preb])
  }
  
  ind.il.lr <- pop.t$curr.state == "il.lr"
  if (sum(ind.il.lr) != 0){
    trans.prob.matrix[ , ind.il.lr]  <- rbind(rep(0, sum(ind.il.lr)), 
                                              1 - p.il.lr2il.hr - p.il.lr2inact - mor.rate[ind.il.lr] - od.rate[ind.il.lr],
                                              rep(p.il.lr2il.hr, sum(ind.il.lr)),
                                              rep(p.il.lr2inact, sum(ind.il.lr)), 
                                              rep(0, sum(ind.il.lr)), 
                                              rep(0, sum(ind.il.lr)), 
                                              mor.rate[ind.il.lr],
                                              od.rate[ind.il.lr])
  }
  
  ind.il.hr <- pop.t$curr.state == "il.hr"
  if (sum(ind.il.hr) != 0){
    trans.prob.matrix[ , ind.il.hr]  <- rbind(rep(0, sum(ind.il.hr)),
                                              rep(p.il.hr2il.lr, sum(ind.il.hr)),
                                              1 - p.il.hr2il.lr - p.il.hr2inact - mor.rate[ind.il.hr] - od.rate[ind.il.hr],
                                              rep(p.il.hr2inact, sum(ind.il.hr)), 
                                              rep(0, sum(ind.il.hr)), 
                                              rep(0, sum(ind.il.hr)), 
                                              mor.rate[ind.il.hr],
                                              od.rate[ind.il.hr])
  }
  
  ind.inact <- pop.t$curr.state == "inact"
  if (sum(ind.inact) != 0){
    trans.prob.matrix[ , ind.inact]  <- rbind(rep(0, sum(ind.inact)),
                                              rep(0, sum(ind.inact)),
                                              rep(0, sum(ind.inact)),
                                              1 - p.inact2relap - mor.rate[ind.inact] - od.rate[ind.inact], 
                                              rep(0, sum(ind.inact)), 
                                              rep(p.inact2relap, sum(ind.inact)), 
                                              mor.rate[ind.inact],
                                              od.rate[ind.inact])
  }
  
  ind.NODU <- pop.t$curr.state == "NODU"
  if (sum(ind.NODU) != 0){
    trans.prob.matrix[ , ind.NODU]   <- rbind(rep(0, sum(ind.NODU)),
                                              rep(0, sum(ind.NODU)),
                                              rep(0, sum(ind.NODU)),
                                              rep(0, sum(ind.NODU)), 
                                              1 - mor.rate[ind.NODU] - od.rate[ind.NODU], 
                                              rep(0, sum(ind.NODU)), 
                                              mor.rate[ind.NODU],
                                              od.rate[ind.NODU])
  }

  ind.relap <- pop.t$curr.state == "relap"
  if (sum(ind.relap) != 0){
    OU.v    <- filter(pop.t, curr.state == "relap")$OU.state
    relap.m <- matrix(0, n.state+1, sum(ind.relap))
    for (r in 1:sum(ind.relap)){
      relap.m[which(OU.v[r] == v.state), r] <- 1 - mor.rate[ind.relap][r] - od.rate[ind.relap][r]
      relap.m[which("dead"  == v.state), r] <- mor.rate[ind.relap][r]          
      relap.m[n.state+1, r]                 <- od.rate[ind.relap][r]
    }
    trans.prob.matrix[ , ind.relap]  <- relap.m
  }
  
  trans.prob.matrix[ , pop.t$curr.state == "dead"]  <- c(0, 0, 0, 0, 0, 0, 1, 0)
  
  return(t(trans.prob.matrix))
}

samplev <- function (probs, m) {
  d <- dim(probs)
  n <- d[1]
  k <- d[2]
  lev <- dimnames(probs)[[2]]
  if (!length(lev)) 
    lev <- 1:k
  ran <- matrix(lev[1], ncol = m, nrow = n)
  U <- t(probs)
  for(i in 2:k) {
    U[i, ] <- U[i, ] + U[i - 1, ]
  }
  if (any((U[k, ] - 1) > 1e-05))
    stop("error in multinom: probabilities do not sum to 1")
  
  for (j in 1:m) {
    un <- rep(runif(n), rep(k, n))
    ran[, j] <- lev[1 + colSums(un > U)]
  }
  ran
}