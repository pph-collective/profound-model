###########################################################################################
#########################       Input and Output setup     ################################
###########################################################################################

###########################################################################################
# This is to define all the input and output variables required in the model (consistent across all scenarios)

# INPUT setup
pop.info    <- c("sex", "race", "age", "residence", "curr.state",
                 "OU.state", "init.age", "init.state", "ever.od", "fx")            # information for each model individual
v.state     <- c("preb", "il.lr", "il.hr", "inact", "NODU", "relap", "dead")       # vector for state names
v.oustate   <- c("preb", "il.lr", "il.hr")                                         # vector for active opioid use state names
n.state     <- length(v.state)                                                     # number of states
n.yr        <- yr.last-yr.first+1
n.t         <- 12 * n.yr                                                           # number of time cycles (in month)
n.rgn       <- length(v.rgn)                                                       # number of regions

# OUTPUT matrices and vectors
v.od        <- rep(0, times = n.t)                                                 # count of overdose events at each time step
v.oddeath   <- rep(0, times = n.t)                                                 # count of overdose deaths at each time step
v.oddeath.w <- rep(0, times = n.t)                                                 # count of overdose deaths that were witnessed at each time step
m.oddeath   <- matrix(0, nrow = n.t, ncol = n.rgn)
colnames(m.oddeath) <- v.rgn
v.odpriv    <- rep(0, times = n.t)                                                 # count of overdose events occurred at private setting at each time step
v.odpubl    <- rep(0, times = n.t)                                                 # count of overdose events occurred at public setting at each time step
v.deathpriv <- rep(0, times = n.t)                                                 # count of overdose deaths occurred at private setting at each time step
v.deathpubl <- rep(0, times = n.t)                                                 # count of overdose deaths occurred at public setting at each time step
v.nlxused   <- rep(0, times = n.t)                                                 # count of naloxone kits used at each time step
v.str       <- c("SQ", "expand", "program")                                        # store the strategy names
cost.item   <- c("TotalCost", "NxCost")
cost.matrix <- matrix(0, nrow=n.t, ncol = length(cost.item))
colnames(cost.matrix) <- cost.item
m.oddeath.fx <- rep(0, times = n.t)                                                # count of overdose deaths with fentanyl present at each time step
m.oddeath.op <- rep(0, times = n.t)                                                # count of overdose deaths among opioid users at each time step
m.oddeath.st <- rep(0, times = n.t)                                                # count of overdose deaths among stimulant users at each time step
m.EDvisits   <- rep(0, times = n.t)                                                # count of opioid overdose-related ED visits at each time step
m.oddeath.hr <- rep(0, times = n.t)                                                # count of overdose deaths among high-risk opioid users (inject heroin) at each time step