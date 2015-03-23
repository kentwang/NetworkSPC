tr = function(A) {
# calculate trace of a square matrix
	return(sum(diag(A)))
}

triadcensus = function(A){
#
# This function computes the triad census as given in Holland and Leinart's
# original paper.  Both directed and undirected networks are allowed.
# 
# Inputs: adjacency matrix A
# Outputs: Triad census proportions T
#
# Todo: Can be refined to remove for loops (Too slow)

	n = dim(A)[1]
	E = matrix(0, n, n)
	M = matrix(0, n, n)
	
	for(i in 1:n) {
		for(j in 1:n) {
			if(A[i, j] == 1 | A[j, i] == 1) {
				E[i, j] = 1
				E[j, i] = 1
			}
			
			if(A[i, j] == 1 && A[j, i] == 1) {
				M[i, j] = 1
				M[j, i] = 1		
			}
		}
	}	

	# Define the following for use in counting configurations
	C=A-M 
	Ec=matrix(1, n, n) - diag(n) - E

	# Counts of 16 configurations in the triad census
	ptm <- proc.time()
	for(k in 1:100){

	T003=tr(Ec %*% Ec %*% Ec)/6;
	T012=sum(sum(Ec %*% Ec * (C + t(C))))/2;
	T102=sum(sum(Ec %*% Ec * M))/2;
	T021D=sum(sum(t(C) %*% C * Ec))/2;
	T021U=sum(sum(C * t(C) * Ec))/2;
	T021C=sum(sum(C %*% C * Ec));
	T111D=sum(sum(A %*% t(A) * Ec - M %*% M * Ec - C %*% t(C) * Ec))/2;
	T111U=sum(sum(t(A) %*% A * Ec - M %*% M * Ec - t(C) %*% C * Ec))/2;
	T030T=sum(sum(C %*% C * C));
	T030C=tr(C %*% C %*% C)/3;
	T201=sum(sum(M %*% M * Ec))/2;
	T120D=sum(sum(t(C) %*% C * M))/2;
	T120U=sum(sum(C %*% t(C) * M))/2;
	T120C=sum(sum(C %*% C * M));
	T210=sum(sum(M %*% M *(C + t(C))))/2;
	T300=tr(M %*% M %*% M)/6;
	}
	proc.time() - ptm


	T=c(T003, T012, T102, T021D, T021U, T021C, T111D, T111U, T030T, T030C, T201, T120D, T120U, T120C, T210, T300)
	
	T = T/choose(dim(A)[1], 3)
	return(T)
}


