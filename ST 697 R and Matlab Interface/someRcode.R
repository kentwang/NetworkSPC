library(ergm)
library(R.matlab)

setwd("../Dropbox/2015 Spring/ST 697/ST 697 R and Matlab Interface/")

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


########### start the simulation for m = 1000  #############
adjR <- list()
for(m in 1:100){
  g.sim = simulate(network(100) ~ transitive + mutual, coef=c(0, 0), nsim = 1)[,]
  cat('iteration: ', m, '\n')
  adjR[[m]] <- g.sim

}
  #adjR = lapply(g.sim, as.matrix.network.adjacency)



  networkSize = length(adjR)
  names(adjR) = paste("mat", 1:networkSize, sep = "")
  
  # write the network to .mat file, read in and reformat into cel array
  writeMat("adjM.mat", adjM = adjR)
  evaluate(matlab, "load('adjM.mat', 'adjM');")
  evaluate(matlab, "adjMSeq = struct2cell(adjM);")

  # temp = getVariable(matlab, "adjM")
  evaluate(matlab, "[mu, Sigma, loadings, theta, uclT2, uclSPE]=phaseIanalysis(adjMSeq, alpha);")

  mu <- getVariable(matlab, "mu")$mu
  Sigma <- getVariable(matlab, "Sigma")$Sigma
  loadings <- getVariable(matlab, "loadings")$loadings
  theta <- getVariable(matlab, "theta")$theta
  uclT2 <- getVariable(matlab, "uclT2")$uclT2
  uclSPE <- getVariable(matlab, "uclSPE")$uclSPE



########### Phase II detection of reciprocity(mutual) = 0.1  #############
N <- 1000
ARL.table <- matrix(NA, 5, 5)
f1 <- 0
for(coef1 in c(0, 0.2, 0.5, 1, 2)){
	f2 <- 0
	f1 <- f1 + 1
	for(coef2 in c(0, 0.2, 0.5, 1, 2)){
	f2 <- f2 + 1

	ARL <- rep(NA, N)
	for(k in 1:N){
	iter <- 0
	repeat{
  	iter <- iter + 1
  	g.sim2 = simulate(network(500) ~ transitive + mutual, coef=c(coef1, coef2), nsim = 1)
  	adjR2 = g.sim2[,]
  	networkSize = length(adjR2)
  	names(adjR2) = paste("mat", 1:networkSize, sep = "")
  
  	# write the network to .mat file, read in
  	writeMat("adjM2.mat", adjM2 = adjR2)
  	evaluate(matlab, "load('adjM2.mat', 'adjM2');")
  	# temp = getVariable(matlab, "adjM2")
 	evaluate(matlab, "[t2, spe]=phaseIImonitoring(adjM2,mu,Sigma,loadings,theta);")
  	t2 <- getVariable(matlab, "t2")$t2
  	spe <- getVariable(matlab, "spe")$spe
  	if((t2 > uclT2) || (spe > uclSPE))break
	}
	ARL[k] <- iter
	}
	ARL.table[f1, f2] <- mean(ARL)
}
}