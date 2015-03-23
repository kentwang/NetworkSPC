library(ergm)
library(network)
library(R.matlab)

#- simulate network data for 1 MC iteration
g.sim = simulate(network(100) ~ transitive + mutual, coef=c(0, 0), nsim = 100)
# adj = lapply(g.sim, FUN = function(v) v[,]) # extract the adjacency matrices as list
adj = lapply(readMat("enron_email_networks_directed.mat")$adj, FUN = function(v) v = v[[1]])

A = adj[[1]]

# Some test code
T = matrix(0, length(adj), 16)
for(i in 1:length(adj)) {
	cat("iter", i, "\n")
	T[i, ] = triadcensus(adj[[i]])	
}

test = lapply(adj, FUN = triadcensus)