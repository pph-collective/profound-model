od.death.target <- c(290, 286, 272, 256)
fx.death.target <- c(0.68, 0.72, 0.83, 0.84)
ed.visit.target <- c(1608, 1679, 1561, 1652)

colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5))

colSums(matrix(sim_sq$m.oddeath.fx, nrow=12, ncol=5)) / colSums(matrix(sim_sq$v.oddeath, nrow=12, ncol=5))

colSums(matrix(sim_sq$m.EDvisits, nrow=12, ncol=5))

colSums(matrix(sim_sq$m.oddeath.st, nrow=12, ncol=5)) / 
  (colSums(matrix(sim_sq$m.oddeath.op, nrow=12, ncol=5)) + colSums(matrix(sim_sq$m.oddeath.st, nrow=12, ncol=5)))

colSums(matrix(sim_sq$m.oddeath.hr, nrow=12, ncol=5))

colSums(matrix(sim_sq$m.oddeath.st, nrow=12, ncol=5))

colSums(matrix(sim_sq$m.oddeath.op, nrow=12, ncol=5))

colSums(matrix(sim_sq$v.deathpriv, nrow=12, ncol=5)) /
  (colSums(matrix(sim_sq$v.deathpriv, nrow=12, ncol=5)) + colSums(matrix(sim_sq$v.deathpubl, nrow=12, ncol=5)))
