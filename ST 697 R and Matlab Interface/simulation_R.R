# This program is a benchmark program
library(ergm)
library(R.matlab)

g.sim.incontrol = simulate(network(100) ~ density + mutual + transitive, 
                           coef=c(0.1, 0.1, 0.1), nsim = 50)
# plot(g.sim.incontrol[[1]])

# justify the network statistics
# summary(g.sim.incontrol[[1]] ~ density + mutual + transitive)
# summary(g.sim.incontrol[[2]] ~ density + mutual + transitive)

g.sim.outcontrol = simulate(network(100) ~ density + mutual + transitive, 
                            coef=c(0.1, 0.1, 1), nsim = 30)

g.sim = c(g.sim.incontrol, g.sim.outcontrol)


# write the sequence into .mat file
adjR = lapply(g.sim, as.matrix.network.adjacency)
networkSize = length(adjR)
names(adjR) = paste("mat", 1:networkSize, sep = "")

writeMat("ergm_sim.mat", adjSeq = adjR)
