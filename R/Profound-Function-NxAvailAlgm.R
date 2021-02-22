nlx.avail  <- function(n.nlx, ou.pop.resid, OD_loc, Low2Priv, nlx.adj){
  n.nlx.loc             <- cbind(n.nlx["high", ] + n.nlx["low", ] * Low2Priv, n.nlx["low", ] * (1-Low2Priv))
  colnames(n.nlx.loc)   <- c("priv", "pub")
  p.nlx.avail           <- n.nlx.loc / (ou.pop.resid$n * t(OD_loc)) * nlx.adj
  return(p.nlx.avail)
}