### This program is an interface calling between R and Matlab for 
### change detection in social network. The interface used is R.matlab in R
### Todo: 
###   - social network simulation using ergm
###      - change the number of vertices as random
###   - unknown reasons for haning on the Matlab server. Code can not run since then.

rm(list = ls())

# install the interface package and ergm
# install.packages("R.matlab")
# install.packages("ergm")
library(R.matlab)
library(ergm)


# start the Matlab server. Needs to pause to get fully booted
Matlab$startServer()
Sys.sleep(5)


# create a Matlab client
matlab = Matlab()

# opent the Matlab client
isOpen = open(matlab)

## The adj matrices should be generated through simulation and passed
# through the interface.
# In the trial, load the adj sequence
# evaluate(matlab, "load('enron_email_networks_directed.mat', 'adj');")

# define the type I error
setVariable(matlab, alpha = 0.05)

# call the existing Matlab functions in the current path. No need to redefine!
# evaluate(matlab, "adj_test = adj(1:10);")
# evaluate(matlab, "[mu, Sigma, loadings, theta, uclT2, uclSPE]=phaseIanalysis(adj_test,alpha);")


## Start the simulation using ergm. 30 in control and 20 out of control
# n.sim = 60
# n.inControl = 40
# n.iter = 0
# for(i in 1:n.sim) {
#   n.iter = n.iter + 1
#   cat(n.iter, "/", n.sim, " iteration\n")
#   
#   # simulate the network in-control or out-of-control
#   if(n.iter <= n.inControl)
#     g.sim = simulate(network(50) ~ edges + mutual, coef=c(0, 0))
#   else
#     g.sim = simulate(network(50) ~ edges + mutual, coef=c(1, 0))
#   
#   # create adjencency matrices sequence in Matlab
#   
# }

# generate network and write into .mat file
#    - looks like we need to save the matrices and combine them in Matlab one by one
# A <- matrix(1:27, ncol=3)
# B <- as.matrix(1:10)
# C <- array(1:18, dim=c(2,3,3))
# D <- list(A, B, A)
# 
# filename="adj_r.mat"
# writeMat(filename, A=A, B=B, C=C, D=D)


########### start the simulation for 100 times #############
n.mc = 100
for(k in 1:n.mc) {
  cat("Monte Carlo iteration", k, "\n")
  
  g.sim = simulate(network(100) ~ edges + mutual, coef=c(0, 0), nsim = 50)
  adjR = lapply(g.sim, as.matrix.network.adjacency)
  networkSize = length(adjR)
  names(adjR) = paste("mat", 1:networkSize, sep = "")
  
  # write the network to .mat file, read in and reformat into cel array
  writeMat("adjM.mat", adjM = adjR)
  evaluate(matlab, "load('adjM.mat', 'adjM');")
  evaluate(matlab, "adjMSeq = struct2cell(adjM);")
  
  # temp = getVariable(matlab, "adjM")
  evaluate(matlab, "[mu, Sigma, loadings, theta, uclT2, uclSPE]=phaseIanalysis(adjMSeq, alpha);")
}



