function Z=linearcombs(A,mu0,Sigma0)
% 
% Function compute the 9 linear combinations of the triad census as given
% in the Holland and Leinart paper.
%
% Inputs:
%   A = adjacency matrix
%   mu0 = 16x1 mean vector of triad census
%   Sigma0 = 16x16 variance-covariance matrix of triad census
%
% Output:
%   Z = 1x9 vector of linear combination of triad census proportions
%
% Weight vectors as given in Holland and Leinart
w=[0,0,3,0,0,0,0,0,0;
   0,1,2,1,0,0,0,0,0;
   1,0,2,2,0,0,0,0,0;
   0,2,1,2,0,1,0,0,0;
   0,2,1,2,1,0,0,0,0;
   0,2,1,2,0,0,0,1,0;
   1,1,1,3,1,0,0,1,0;
   1,1,1,3,0,1,0,1,1;
   0,3,0,3,1,1,1,0,0;
   0,3,0,3,0,0,0,3,0;
   2,0,1,4,1,1,0,2,2;
   1,2,0,4,2,1,2,0,0;
   1,2,0,4,1,2,2,0,0;
   1,2,0,4,1,1,1,2,1;
   2,1,0,5,2,2,3,1,1;
   3,0,0,6,3,3,6,0,0];
%
% Compute triad census (proportions)
T=triadcensus(A)';
%
% Compute linear combinations of triad census proportions
Z2=(w(:,1)'*T-w(:,1)'*mu0)/sqrt(w(:,1)'*Sigma0*w(:,1));
Z3=(w(:,2)'*T-w(:,2)'*mu0)/sqrt(w(:,2)'*Sigma0*w(:,2));
Z4=(w(:,3)'*T-w(:,3)'*mu0)/sqrt(w(:,3)'*Sigma0*w(:,3));
Z5=(w(:,4)'*T-w(:,4)'*mu0)/sqrt(w(:,4)'*Sigma0*w(:,4));
Z6=(w(:,5)'*T-w(:,5)'*mu0)/sqrt(w(:,5)'*Sigma0*w(:,5));
Z7=(w(:,6)'*T-w(:,6)'*mu0)/sqrt(w(:,6)'*Sigma0*w(:,6));
Z8=(w(:,7)'*T-w(:,7)'*mu0)/sqrt(w(:,7)'*Sigma0*w(:,7));
Z9=(w(:,8)'*T-w(:,8)'*mu0)/sqrt(w(:,8)'*Sigma0*w(:,8));
Z10=(w(:,9)'*T-w(:,9)'*mu0)/sqrt(w(:,9)'*Sigma0*w(:,9));
%
% Store in 1 x 9 vector
Z=[Z2 Z3 Z4 Z5 Z6 Z7 Z8 Z9 Z10];
%
%
