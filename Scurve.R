x <- c(2012:2025)

b <- 2020

inv_logit <- function(x) {
  return(0.5 / (1 + exp(- (x - b))))
}

y <- inv_logit(x)

plot(y ~ x)