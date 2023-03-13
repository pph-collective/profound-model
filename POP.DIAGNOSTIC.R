ss=3
seed <- seed.sameple[ss]
dead_pop <- NULL

sim_sq <- MicroSim(init_ppl, params = params, timesteps, agent_states, discount.rate, PT.out = T, strategy = "SQ", seed = seed)        # run for status quo
pop.trace <- sim_sq$pop.trace
length(pop.trace)
fun <- function(DT){
  return(data.frame(table(DT$curr.state)))
}

agent_states <- agent_states[agent_states !="dead"]
pop.dist <- lapply(pop.trace, FUN = fun)
pop.dist.tb <- matrix(0, ncol = length(agent_states), nrow = length(pop.dist))
rownames(pop.dist.tb) <- c(1:length(pop.dist))
colnames(pop.dist.tb) <- agent_states
pop.long.table <- data.frame(month = rep(1:length(pop.trace), each = 3*length(agent_states)), race = rep(rep(c("white", "black", "hisp"), each = length(agent_states)), length(pop.trace)), state = rep(agent_states, 3 * length(pop.trace)), value = 0)
pop.ever.od.table <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
pop.fx.table <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
ever.odamongil.hr <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
fxamongil.hr      <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
il.hr.no          <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
il.hramongtottal  <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
ever.odamongil.hr.abs <- data.frame(month = rep(1:length(pop.trace), each = 3), race = rep(c("white", "black", "hisp"), length(pop.trace)), value = 0)
for(i in 1:length(pop.dist)){
  pop.dist.tb[i, "preb"] <- pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "preb"]
  pop.dist.tb[i, "il.lr"] <- pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "il.lr"]
  pop.dist.tb[i, "il.hr"] <- pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "il.hr"]
  pop.dist.tb[i, "inact"] <- pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "inact"]
  pop.dist.tb[i, "NODU"] <- pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "NODU"]
  pop.dist.tb[i, "relap"] <- pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "relap"]
  # pop.dist.tb[i, "dead"] <- ifelse(sum(pop.dist[[i]]$Var1 == "dead") == 1, pop.dist[[i]]$Freq[pop.dist[[i]]$Var1 == "dead"],0)
  table.state <- table(pop.trace[[i]]$race, pop.trace[[i]]$curr.state)
  table.state <- table.state[ , agent_states]
  pop.long.table$value[pop.long.table$month == i & pop.long.table$race == "white"] <- table.state["white",]
  pop.long.table$value[pop.long.table$month == i & pop.long.table$race == "black"] <- table.state["black",]
  pop.long.table$value[pop.long.table$month == i & pop.long.table$race == "hisp"] <- table.state["hisp",]
  pop.ever.od.table$value[pop.ever.od.table$month == i & pop.ever.od.table$race == "white"] <- sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state != "dead"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state != "dead", ])
  pop.ever.od.table$value[pop.ever.od.table$month == i & pop.ever.od.table$race == "black"] <- sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state != "dead"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state != "dead", ])
  pop.ever.od.table$value[pop.ever.od.table$month == i & pop.ever.od.table$race == "hisp"] <- sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state != "dead"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state != "dead", ])
  pop.fx.table$value[pop.fx.table$month == i & pop.fx.table$race == "white"] <- sum(pop.trace[[i]]$fx[pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state != "dead"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state != "dead", ])
  pop.fx.table$value[pop.fx.table$month == i & pop.fx.table$race == "black"] <- sum(pop.trace[[i]]$fx[pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state != "dead"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state != "dead", ])
  pop.fx.table$value[pop.fx.table$month == i & pop.fx.table$race == "hisp"]  <- sum(pop.trace[[i]]$fx[pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state != "dead"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state != "dead", ])
  ever.odamongil.hr$value[ever.odamongil.hr$month == i & ever.odamongil.hr$race == "white"] <-  sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr", ])
  ever.odamongil.hr$value[ever.odamongil.hr$month == i & ever.odamongil.hr$race == "black"] <-  sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr", ])
  ever.odamongil.hr$value[ever.odamongil.hr$month == i & ever.odamongil.hr$race == "hisp"]  <-  sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr", ])
  ever.odamongil.hr.abs$value[ever.odamongil.hr.abs$month == i & ever.odamongil.hr.abs$race == "white"] <-  sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr"])
  ever.odamongil.hr.abs$value[ever.odamongil.hr.abs$month == i & ever.odamongil.hr.abs$race == "black"] <-  sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr"])
  ever.odamongil.hr.abs$value[ever.odamongil.hr.abs$month == i & ever.odamongil.hr.abs$race == "hisp"]  <-  sum(pop.trace[[i]]$ever.od[pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr"])
  fxamongil.hr$value[fxamongil.hr$month == i & fxamongil.hr$race == "white"] <-  sum(pop.trace[[i]]$fx[pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr", ])
  fxamongil.hr$value[fxamongil.hr$month == i & fxamongil.hr$race == "black"] <-  sum(pop.trace[[i]]$fx[pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr", ])
  fxamongil.hr$value[fxamongil.hr$month == i & fxamongil.hr$race == "hisp"]  <-  sum(pop.trace[[i]]$fx[pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr"]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr", ])
  il.hr.no$value[il.hr.no$month == i & il.hr.no$race == "white"] <-  nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr", ])
  il.hr.no$value[il.hr.no$month == i & il.hr.no$race == "black"] <-  nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr", ])
  il.hr.no$value[il.hr.no$month == i & il.hr.no$race == "hisp"]  <-  nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr", ])
  il.hramongtottal$value[il.hramongtottal$month == i & il.hramongtottal$race == "white"] <-  nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state == "il.hr", ]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "white" & pop.trace[[i]]$curr.state != "dead", ])
  il.hramongtottal$value[il.hramongtottal$month == i & il.hramongtottal$race == "black"] <-  nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state == "il.hr", ]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "black" & pop.trace[[i]]$curr.state != "dead", ])
  il.hramongtottal$value[il.hramongtottal$month == i & il.hramongtottal$race == "hisp"]  <-  nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state == "il.hr", ]) / nrow(pop.trace[[i]][pop.trace[[i]]$race == "hisp" & pop.trace[[i]]$curr.state != "dead", ])
} 

calib.comp$white_sim[calib.comp$par == "ODdeaths"] <- colSums(matrix(rowSums(sim_sq$m.oddeath.white), nrow = 12))[1:5]
calib.comp$black_sim[calib.comp$par == "ODdeaths"] <- colSums(matrix(rowSums(sim_sq$m.oddeath.black), nrow = 12))[1:5]
calib.comp$hisp_sim[calib.comp$par == "ODdeaths"]  <- colSums(matrix(rowSums(sim_sq$m.oddeath.hisp), nrow = 12))[1:5]
calib.comp$white_sim[calib.comp$par == "Fx_ODD"] <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "white"], nrow = 12))[1:5] / calib.comp$white_sim[calib.comp$par == "ODdeaths"]
calib.comp$black_sim[calib.comp$par == "Fx_ODD"] <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "black"], nrow = 12))[1:5] / calib.comp$black_sim[calib.comp$par == "ODdeaths"]
calib.comp$hisp_sim[calib.comp$par == "Fx_ODD"]  <- colSums(matrix(sim_sq$m.oddeath.fx.race[, "hisp"], nrow = 12))[1:5] / calib.comp$hisp_sim[calib.comp$par == "ODdeaths"]
calib.comp$white_sim[calib.comp$par == "EDvisits"] <- colSums(matrix(sim_sq$m.EDvisits.race[, "white"], nrow = 12))[1:5]
calib.comp$black_sim[calib.comp$par == "EDvisits"] <- colSums(matrix(sim_sq$m.EDvisits.race[, "black"], nrow = 12))[1:5]
calib.comp$hisp_sim[calib.comp$par == "EDvisits"]  <- colSums(matrix(sim_sq$m.EDvisits.race[, "hisp"], nrow = 12))[1:5]

calib.comp
colSums(matrix(rowSums(sim_sq$m.oddeath), nrow = 12))
pop.dist.tb

library(ggplot2)
pop.dist.tf <- reshape2::melt(pop.dist.tb, value.name = "Number")
names(pop.dist.tf) <- c("Month", "State", "Number")

ggplot(data=pop.dist.tf, aes(x=Month, y=Number, group=State)) +
  geom_line(aes(col = State))+
  geom_point(aes(col = State)) +
  ggtitle("State distribution")

ggplot(data=pop.long.table, aes(x=month, y=value, group=state)) +
  geom_line(aes(col = state))+
  geom_point(aes(col = state)) +
  facet_wrap( ~ race, nrow = 1, scales = "free") +
  ggtitle("State distribution by race")

ggplot(data=pop.ever.od.table, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race)) + 
  ggtitle("ever od by race")

ggplot(data=pop.fx.table, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race)) +
  ggtitle("fentanyl exposure by race")

ggplot(data=ever.odamongil.hr, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race)) +
  ggtitle("ever od among illicit high-risk by race")

ggplot(data=ever.odamongil.hr.abs, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race)) +
  ggtitle("ever od among illicit high-risk absolute number by race")

ggplot(data=fxamongil.hr, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race))+
  ggtitle("fentanyl exposure among illicit high-risk by race")

ggplot(data=il.hr.no, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race)) +
  ggtitle("illicit high-risk absolute number by race")

ggplot(data=il.hramongtottal, aes(x=month, y=value, group=race)) +
  geom_line(aes(col = race))+
  geom_point(aes(col = race)) +
  ggtitle("illicit high-risk proportion by race")


# dead.pop<- pop.trace[[96]][pop.trace[[96]]$curr.state == "dead",]
dead_pop <- sim_sq$dead_pop
table(dead_pop$OU.state, dead_pop$fx)
table(dead_pop$OU.state, dead_pop$ever.od)
table(dead_pop$race, dead_pop$tmi)/rep(colSums(table(dead_pop$race, dead_pop$tmi)), each =3)
table(dead_pop$OU.state, dead_pop$tmi)
table(dead_pop$fx, dead_pop$tmi)
table(dead_pop$ind)


dead.ind<-as.numeric(row.names(as.matrix(table(dead_pop$ind))))
dead.mt<-which(table(dead_pop$ind)>6)
table(dead_pop[dead_pop$ind %in% dead.ind[dead.mt],]$residence)
dead_pop[dead_pop$ind %in% dead.ind[which(table(dead_pop$ind)==max(table(dead_pop$ind)))],]