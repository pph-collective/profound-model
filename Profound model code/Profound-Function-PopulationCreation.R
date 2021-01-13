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
#################        Population creation function      #############################
########################################################################################

pop.creat <- function(pop, mean.age, pop.prop.list, seed = 2020){
  list2env(pop.prop.list, environment())
  for (i in 1:nrow(pop)){
    set.seed(seed+i)
    # determine age
    age        <- init.age <- rpois(1, mean.age)
    
    # determine residence
    rds        <- colnames(prop.rsd)
    residence  <- rds[sample(1:length(rds), size = 1, prob = prop.rsd[1,])]
    
    # determine drug use state
    du.state   <- v.state[1:5]
    du.prob    <- c(prop.preb, prop.illicit * prop.il.lr, prop.illicit * prop.il.hr, prop.inact, prop.NOUD)
    curr.state <- init.state <- du.state[sample(1:length(du.state), size = 1, prob = du.prob)]
    if (curr.state == "inact"){
      OU.state <- du.state[sample(1:length(du.state[1:3]), size = 1, prob = du.prob[1:3]/sum(du.prob[1:3]))]
    } else {
      OU.state <- curr.state
    }
    
    # determine fentanyl use
    fx         <- sample(0:1, size = 1, prob = c(1-prop.fx, prop.fx))
    
    pop$age[i]        <- age
    pop$residence[i]  <- residence
    pop$curr.state[i] <- curr.state
    pop$OU.state[i]   <- OU.state
    pop$init.age[i]   <- init.age
    pop$init.state[i] <- init.state
    pop$fx[i]         <- fx
    pop$ever.od[i]    <- 0
  }
  return(pop)
}