# TODO file header

nlx.avail.algm  <- function(n.nlx, ou.pop.resid, OD_loc, Low2Priv, nlx.adj, cap){
  n.nlx.loc             <- cbind(n.nlx["high", ] + n.nlx["low", ] * Low2Priv, n.nlx["low", ] * (1-Low2Priv))
  colnames(n.nlx.loc)   <- c("priv", "pub")
  # p.nlx.avail           <- n.nlx.loc / (ou.pop.resid$n * t(OD_loc)) * nlx.adj
  p.nlx.avail           <- cap * (1 - exp(-(nlx.adj / cap) * n.nlx.loc / (ou.pop.resid$n * t(OD_loc))))
  return(p.nlx.avail)
}