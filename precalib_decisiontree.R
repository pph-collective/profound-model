###############################################################################################
###################### PROFOUND Naloxone Distribution model #### 2020 #########################
###############################################################################################
# A pre-step for calibrating decision tree module to ensure fitting to the following targets:
# 1. % of non-fatal opioid overdoses resulting in EMS (911 called) in public settings - 30.8%, assumed 95% CI: [24.6%, 37.0%] -> se: 0.0314
# 2. % of fatal opioid overdoses occurring in public settings - 10%, assumed 95% CI: [8%-12%] -> se:0.0102  
#
# Through adjusting the following model parameters:
# 1. Proportion of overdoses occurring in public versus private settings
# 2. Multiplier/ratio of probability of overdose being witnessed in private versus public settings
# 3. Multiplier/ratio of probability of witnessed overdose resulting in 911 called in private versus public settings
#
# All other parameters using point estimates
# The calibration is performed using a Bayesian method
# We aim to derived 1,000,000 acceptable sets for subsequent analysis 
#
###############################################################################################
#########################           Decision Tree           ###################################
###############################################################################################
#### This decision tree module has been modified to disable random sampling at forks       ####
#### Instead of individual-based, the modified decision tree is cohort based               ####
#### Many parameters have been simplified - yet the primary results should remain similar  ####
###############################################################################################
#############################################################################
# 1. Decision tree parameters
#############################################################################
# INPUT PARAMETERS and Packages

rm(list = ls())
library("argparser")
library(openxlsx)
library(abind)
# calibration functionality
library(lhs)
library(IMIS)
library(matrixStats) # package used for sumamry statistics
# visualization
library(plotrix)
library(psych)
source("data_input.R")
fix.params <- list()
fix.params$OD_wit_pub  <- params$OD_wit_pub
fix.params$OD_911_pub  <- params$OD_911_pub
fix.params$OD_hosp     <- params$OD_hosp
fix.params$mor_bl      <- params$mor_bl
fix.params$rr_mor_EMS  <- params$rr_mor_EMS
fix.params$mor_nx      <- params$mor_nx

## Customizable function inputs 
n.od          <- 10000  # number of overdoses
nlx.avai.pub  <- 0.1    # naloxone available in public settings
nlx.avai.priv <- 0.15   # naloxone available in private settings

#############################################################################
# 2. Decision tree function
#############################################################################

decision_tree <- function(n.od, nlx.avai.pub, nlx.avai.priv, fix.params, cal.params) {
  list2env(fix.params, environment())
  # out.colnames <- c("od.death", "EMS", "locpriv", "nlx.used", "wtns")
  # decntree.out <- matrix(0, nrow = n.od, ncol = length(out.colnames))
  # colnames(decntree.out) <- out.colnames
  
  OD_pub <- cal.params["OD_pub"]
  rr_OD_wit_priv <- cal.params["rr_OD_wit_priv"]
  rr_OD_911_priv <- cal.params["rr_OD_911_priv"]
  
  pp1 <- OD_pub * OD_wit_pub * nlx.avai.pub * OD_911_pub * OD_hosp * (1-mor_nx)
  pp2 <- OD_pub * OD_wit_pub * nlx.avai.pub * OD_911_pub * OD_hosp * mor_nx
  
  pp3 <- OD_pub * OD_wit_pub * nlx.avai.pub * OD_911_pub * (1-OD_hosp) * (1-mor_nx) 
  pp4 <- OD_pub * OD_wit_pub * nlx.avai.pub * OD_911_pub * (1-OD_hosp) * mor_nx
  
  pp5 <- OD_pub * OD_wit_pub * nlx.avai.pub * (1- OD_911_pub) * (1-mor_nx)
  pp6 <- OD_pub * OD_wit_pub * nlx.avai.pub * (1- OD_911_pub) * mor_nx
  
  pp7 <- OD_pub * OD_wit_pub * (1-nlx.avai.pub) * OD_911_pub * OD_hosp * (1-mor_bl*rr_mor_EMS)
  pp8 <- OD_pub * OD_wit_pub * (1-nlx.avai.pub) * OD_911_pub * OD_hosp * (mor_bl*rr_mor_EMS)
  
  pp9 <- OD_pub * OD_wit_pub * (1-nlx.avai.pub) * OD_911_pub * (1-OD_hosp) * (1-mor_bl*rr_mor_EMS)
  pp10<- OD_pub * OD_wit_pub * (1-nlx.avai.pub) * OD_911_pub * (1-OD_hosp) * (mor_bl*rr_mor_EMS)
  
  pp11<- OD_pub * OD_wit_pub * (1-nlx.avai.pub) * (1- OD_911_pub) * (1-mor_bl)
  pp12<- OD_pub * OD_wit_pub * (1-nlx.avai.pub) * (1- OD_911_pub) * mor_bl
  
  pp13<- OD_pub * (1- OD_wit_pub) * (1-mor_bl)
  pp14<- OD_pub * (1- OD_wit_pub) * mor_bl
  
  pp15<- (1-OD_pub) * (1- OD_wit_pub * rr_OD_wit_priv) * (1-mor_bl)
  pp16<- (1-OD_pub) * (1- OD_wit_pub * rr_OD_wit_priv) * mor_bl
  
  pp17<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * nlx.avai.priv * (OD_911_pub * rr_OD_911_priv) * OD_hosp * (1-mor_nx)
  pp18<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * nlx.avai.priv * (OD_911_pub * rr_OD_911_priv) * OD_hosp * mor_nx
  
  pp19<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * nlx.avai.priv * (OD_911_pub * rr_OD_911_priv) * (1-OD_hosp) * (1-mor_nx)
  pp20<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * nlx.avai.priv * (OD_911_pub * rr_OD_911_priv) * (1-OD_hosp) * mor_nx
  
  pp21<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * nlx.avai.priv * (1-OD_911_pub * rr_OD_911_priv) * (1-mor_nx)
  pp22<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * nlx.avai.priv * (1-OD_911_pub * rr_OD_911_priv) * mor_nx
  
  pp23<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * (1-nlx.avai.priv) * (OD_911_pub * rr_OD_911_priv) * OD_hosp * (1-mor_bl*rr_mor_EMS)
  pp24<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * (1-nlx.avai.priv) * (OD_911_pub * rr_OD_911_priv) * OD_hosp * (mor_bl*rr_mor_EMS)
  
  pp25<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * (1-nlx.avai.priv) * (OD_911_pub * rr_OD_911_priv) * (1-OD_hosp) * (1-mor_bl*rr_mor_EMS)
  pp26<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * (1-nlx.avai.priv) * (OD_911_pub * rr_OD_911_priv) * (1-OD_hosp) * (mor_bl*rr_mor_EMS)
  
  pp27<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * (1-nlx.avai.priv) * (1-OD_911_pub * rr_OD_911_priv) * (1-mor_bl)
  pp28<- (1-OD_pub) * (OD_wit_pub * rr_OD_wit_priv) * (1-nlx.avai.priv) * (1-OD_911_pub * rr_OD_911_priv) * mor_bl

  np <- n.od * c(pp1, pp2, pp3, pp4, pp5, pp6, pp7, pp8, pp9, pp10, pp11, pp12, pp13, pp14, pp15, pp16, pp17, pp18, pp19, pp20, pp21, pp22, pp23, pp24, pp25, pp26, pp27, pp28)
  
  nonfatal_ems_pub   <- sum(np[c(1, 3, 7, 9)])
  nonfatal_ems_priv  <- sum(np[c(17, 19, 23, 25)])
  per_nonfatal_ems_pub <- nonfatal_ems_pub /(nonfatal_ems_pub + nonfatal_ems_priv)
  
  fatal_pub  <- sum(np[c(2, 4, 6, 8, 10, 12, 14)])
  fatal_priv <- sum(np[c(16, 18, 20, 22, 24, 26, 28)])
  per_fatal_pub <- fatal_pub /(fatal_pub + fatal_priv)
  
  return(c(per_nonfatal_ems_pub, per_fatal_pub))
}

#### Sample priors ####
# names and number of input parameters to be calibrated
param.names <- c("OD_pub", "rr_OD_wit_priv", "rr_OD_911_priv")
n.param <- length(param.names)
# range on input search space
lb <- c(OD_pub = 0.05, rr_OD_wit_priv = 0.2, rr_OD_911_priv = 0.4) # lower bound
ub <- c(OD_pub = 0.31, rr_OD_wit_priv = 1, rr_OD_911_priv = 0.7) # upper bound

sample.prior <- function(n.samp){
  m.lhs.unit   <- lhs::randomLHS(n = n.samp, k = n.param)
  m.param.samp <- matrix(nrow = n.samp, ncol = n.param)
  colnames(m.param.samp) <- param.names
  for (i in 1:n.param){
    m.param.samp[, i] <- qunif(m.lhs.unit[,i],
                               min = lb[i],
                               max = ub[i])
  }
  return(m.param.samp)
}

# Prior density
f_log_prior <- function(v.params){
  if(is.null(dim(v.params))) { # If vector, change to matrix
    v.params <- t(v.params) 
  }
  n.samp <- nrow(v.params)
  colnames(v.params) <- param.names
  lprior <- rep(0, n.samp)
  for (i in 1:n.param){
    lprior <- lprior + dunif(v.params[, i],
                             min = lb[i],
                             max = ub[i], 
                             log = T)
  }
  return(lprior)
}

prior <- function(v.params) { 
  exp(f_log_prior(v.params)) 
}

# Likelihood function
# number of calibration targets
## Additional Targets
n.target <- 2  # or any number of targets that you have

# define log-likelihood function
f_llik_2tar <- function(v.params){
  # par_vector: a vector (or matrix) of model parameters 
  if(is.null(dim(v.params))) { # If vector, change to matrix
    v.params <- t(v.params) 
  }
  n.samp <- nrow(v.params)
  v.llik <- matrix(0, nrow = n.samp, ncol = n.target) 
  llik.overall <- numeric(n.samp)
  for(j in 1:n.samp) { # j=1
    jj <- tryCatch( { 
      ###   Run model for parametr set "v.params" ###
      model.res <- decision_tree(n.od= n.od, nlx.avai.pub = nlx.avai.pub, nlx.avai.priv = nlx.avai.priv, fix.params = fix.params, v.params[j, ])
      
      ###  Calculate log-likelihood of model outputs to targets  ###
      # TARGET 1: perc_nonfatal_ems_pub
      # log likelihood  
      v.llik[j, 1] <- sum(dnorm(x = 0.308,
                                mean = model.res[1],
                                sd = 0.0314,
                                log = T))
      
      # TARGET 2: perc_fatal_pub
      # log likelihood
      v.llik[j, 2] <- sum(dnorm(x = 0.1,
                                mean = model.res[2],
                                sd = 0.0102,
                                log = T))
      
      # OVERALL 
      llik.overall[j] <- sum(v.llik[j, ])
    }, error = function(e) NA) 
    if(is.na(jj)) { llik.overall <- -Inf }
  } # End loop over sampled parameter sets
  # return LLIK
  return(llik.overall)
}

# define likelihood function
likelihood <- function(v.params){ 
  exp(f_llik_2tar(v.params)) 
}

# Log-posterior
f_log_post <- function(v.params) { 
  lpost <- f_log_prior(v.params) + f_llik(v.params)
  return(lpost) 
}

# calibration
# number of resamples
n.resamp <- 1000000

# run IMIS
fit.imis <- IMIS(B = 1000, # the incremental sample size at each iteration of IMIS.
                 B.re = n.resamp, # the desired posterior sample size
                 number_k = 100, # the maximum number of iterations in IMIS.
                 D = 0)
# obtain posterior
m.calib.post <- fit.imis$resample
head(m.calib.post)

# # Plot the 1000 draws from the posterior with marginal histograms
# psych::pairs.panels(m.calib.post)

# Compute posterior mean
v.calib.post.mean <- colMeans(m.calib.post)
v.calib.post.mean

# Compute posterior median and 95% credible interval
m.calib.post.95cr <- matrixStats::colQuantiles(m.calib.post, probs = c(0.025, 0.5, 0.975))
m.calib.post.95cr

# compute maximum a posteriori
v.calib.like <- likelihood(m.calib.post)
v.calib.post.map <- m.calib.post[which.max(v.calib.like), ]

# Calibrated model predictions vs. targets
v.out.post.map <- decision_tree(n.od= n.od, nlx.avai.pub = nlx.avai.pub, nlx.avai.priv = nlx.avai.priv, fix.params = fix.params, v.calib.post.map)
perc_nonfatal_ems_pub <- perc_fatal_pub <- numeric(n.resamp)
for (i in 1:n.resamp){
  rests <- decision_tree(n.od= n.od, nlx.avai.pub = nlx.avai.pub, nlx.avai.priv = nlx.avai.priv, fix.params = fix.params, m.calib.post[i,])
  perc_nonfatal_ems_pub[i] <- rests[1]
  perc_fatal_pub[i] <- rests[2]
}

library(gridExtra)
library(ggplot2)
library(reshape2)
post.params <- data.frame(OD_pub = m.calib.post[, "OD_pub"], rr_OD_wit_priv = m.calib.post[, "rr_OD_wit_priv"], rr_OD_911_priv = m.calib.post[, "rr_OD_911_priv"])
post.params$sample <- c(1:nrow(m.calib.post))
post.params <- melt(post.params, id.vars = "sample")
post.params$variable <- factor(post.params$variable, levels = c("OD_pub", "rr_OD_wit_priv", "rr_OD_911_priv"), 
                  labels = c("Proportion of overdoses occurring in public", "RR of overdose witnessed in private/semi-private vs public", "RR of witness calling EMS in private/semi-private vs public"))
fig1 <- ggplot(post.params, aes(x=value)) + geom_histogram() + labs(y= "Count number", x = "Posterior distribution of parameters") + facet_wrap(~variable, nrow = 1, scales = "free")
fig2 <- ggplot(as.data.frame(perc_nonfatal_ems_pub, aes), aes(x = perc_nonfatal_ems_pub)) + geom_histogram(color="black", fill="white") + labs(y= "Count number", x = "Percentage of non-fatal opioid overdoses resulting in EMS in public settings") +
  geom_vline(xintercept = 0.308, colour = "red", size =2)
fig3 <- ggplot(as.data.frame(perc_fatal_pub, aes), aes(x = perc_fatal_pub)) + geom_histogram(color="black", fill="white") + labs(y= "Count number", x = "Percentage of fatal opioid overdoses in public settings") +
  geom_vline(xintercept = 0.1, colour = "red", size =2)

grid.arrange(
  fig1, fig2, fig3,
  widths = c(1, 1, 1, 1, 1, 1),
  layout_matrix = rbind(c(1, 1, 1, 1, 1, 1),
                        c(2, 2, 2, 3, 3, 3))
)

saveRDS(m.calib.post, file = "Inputs/Precalib_dtree.rds")
