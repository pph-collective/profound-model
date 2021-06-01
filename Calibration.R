#Proportion of OD deaths atributable to opioid use
colSums(matrix(sim_sq$m.oddeath.op, nrow = 12)) / (colSums(matrix(sim_sq$m.oddeath.op, nrow = 12)) + colSums(matrix(sim_sq$m.oddeath.st, nrow = 12)))

#Proportion of OD deaths atributable to stimulant use
colSums(matrix(sim_sq$m.oddeath.st, nrow = 12)) / (colSums(matrix(sim_sq$m.oddeath.op, nrow = 12)) + colSums(matrix(sim_sq$m.oddeath.st, nrow = 12)))

#Total number of opioid OD deaths at the state level
(colSums(matrix(sim_sq$m.oddeath.op, nrow = 12)) + colSums(matrix(sim_sq$m.oddeath.st, nrow = 12)))

#Proportion of OD deaths with fentanyl present
colSums(matrix(sim_sq$m.oddeath.fx, nrow = 12)) / (colSums(matrix(sim_sq$m.oddeath.op, nrow = 12)) + colSums(matrix(sim_sq$m.oddeath.st, nrow = 12)))

#ED visits
colSums(matrix(sim_sq$m.EDvisits, nrow = 12))

#OD deaths at the city/town level in 2019
data <- data.frame(City = colnames(sim_sq$m.oddeath),
  ODdeaths =colSums(sim_sq$m.oddeath[ 37:48, ]))
library(ggplot2)
p<-ggplot(data=data, aes(x=reorder(City, desc(City)), y=ODdeaths)) +
  xlab("City/Town") + ylab("Opioid OD deaths") +
  geom_bar(stat="identity") +
  coord_flip()
p