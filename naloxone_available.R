# TODO file header

nlx.avail.algm <- function(nlx_avail, ou.pop.resid, OD_loc, Low2Priv, nlx.adj, cap) {
  nlx_avail_loc <- cbind(nlx_avail["high", ] + nlx_avail["low", ] * Low2Priv, nlx_avail["low", ] * (1 - Low2Priv))
  colnames(nlx_avail_loc) <- c("priv", "pub")
  p.nlx.avail <- cap * (1 - exp(-(nlx.adj / cap) * nlx_avail_loc / (ou.pop.resid$n * t(OD_loc))))
  return(p.nlx.avail)
}
