
# Histogram showing variability in individual total costs
plot(density(tc), main = paste("Total cost per person"), xlab = "Cost ($)")

# Histogram showing variability in individual total QALYs
plot(density(te), main = paste("Total QALYs per person"), xlab = "QALYs")


# PLOT THE DISTRIBUTION OF THE POPULATON ACROSS HEALTH STATES OVER TIME (TRACE)
# count the number of individuals in each health state at each cycle
m.TR <- t(apply(m.M, 2, function(x) table(factor(x, levels = v.n, ordered = TRUE)))) 
m.TR <- m.TR / n.i                                       # calculate the proportion of individuals 
colnames(m.TR) <- v.n                                    # name the rows of the matrix
rownames(m.TR) <- paste("Cycle", 0:n.t, sep = " ")       # name the columns of the matrix

# Plot trace of first health state
plot(0:n.t, m.TR[, 1], type = "l", main = "Health state trace", 
     ylim = c(0, 1), ylab = "Proportion of cohort", xlab = "Cycle")
# add a line for each additional state
for (n.s in 2:length(v.n)) {
  lines(m.TR[, n.s], col = n.s)       # adds a line to current plot
}
legend("topright", v.n, col = 1:3,    # add a legend to current plot
       lty = rep(1, 3), bty = "n", cex = 0.65)


