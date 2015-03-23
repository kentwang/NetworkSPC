library(ergm)
library(R.matlab)

# setwd("../Documents/GitHub/NetworkSPC/ST 697 R and Matlab Interface/")

# test the function
# test = NetworkSPC_MC(150, 100, 600, c(0, 0.2, 0.5, 1, 2), c(0, 0.2, 0.5, 1, 2), 0.05)

NetworkSPC_MC = function(M, N, meanNumNodes, transitivity, mutuality, alpha) {
	
	# initiate the list to be returned. 
	#- todo: May include network information
	result = list()	

	currentTime = proc.time()

	# start the Matlab server. Needs to pause to get fully booted
	Matlab$startServer()
	Sys.sleep(5)


	# create a Matlab client. Don't need to create every time
	matlab = Matlab()

	# opent the Matlab client
	isOpen = open(matlab)
	if(isOpen) cat("Matlab client is connected!\n")
	
	##- setting related parameters	
	# define the type I error
	setVariable(matlab, alpha = 0.05)
	
	# generate number of nodes. Using Uniform for now. 
	#- todo: May use a sequence as input (Enron)
	numNodes = round(runif(M, 200, 2 * meanNumNodes - 200))

	##- start the simulation for M -##
	adjR <- list()
	for(m in 1:M){
		g.sim = simulate(network(numNodes[m]) ~ transitive + mutual, coef=c(0, 0), nsim = 1)[,]
		cat('iteration: ', m, '\n')
		adjR[[m]] <- g.sim
	}
		
	networkSize = length(adjR)
	names(adjR) = paste("mat", 1:networkSize, sep = "")
	
	# write the network to .mat file, read in and reformat into cel array
	writeMat("adjM.mat", adjM = adjR)
	evaluate(matlab, "load('adjM.mat', 'adjM');")
	evaluate(matlab, "adjMSeq = struct2cell(adjM);")

	evaluate(matlab, "[mu, Sigma, loadings, theta, uclT2, uclSPE]=phaseIanalysis(adjMSeq, alpha);")

	mu <- getVariable(matlab, "mu")$mu
	Sigma <- getVariable(matlab, "Sigma")$Sigma
	loadings <- getVariable(matlab, "loadings")$loadings
	theta <- getVariable(matlab, "theta")$theta
	uclT2 <- getVariable(matlab, "uclT2")$uclT2
	uclSPE <- getVariable(matlab, "uclSPE")$uclSPE



	##- Phase II detection -##
	ARL.table <- matrix(NA, length(transitivity), length(mutuality))
	f1 <- 0
	for(coef1 in transitivity){
		f2 <- 0
		f1 <- f1 + 1
		for(coef2 in mutuality){
			f2 <- f2 + 1

			ARL <- rep(NA, N)
			for(k in 1:N){	
				cat("Transitivity", coef1, "Mutuality", coef2, 'iter', k, "\n")

				iter <- 0
				repeat{
		  			iter <- iter + 1
				  	g.sim2 = simulate(network(100) ~ transitive + mutual, coef=c(coef1, coef2), nsim = 1)
				  	adjR2 = g.sim2[,]
				  	networkSize = length(adjR2)
				  	names(adjR2) = paste("mat", 1:networkSize, sep = "")
  
				  	# write the network to .mat file, read in
				  	writeMat("adjM2.mat", adjM2 = adjR2)
				  	evaluate(matlab, "load('adjM2.mat', 'adjM2');")
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

	runningTime = proc.time() - currentTime

	# close the matlab client
	close(matlab)
	
	result$ARL.table = ARL.table
	result$runningTime = runningTime
	result$numNodes = numNodes
	result$transitivity = transitivity
	result$mutuality = mutuality

	return(result)
}
