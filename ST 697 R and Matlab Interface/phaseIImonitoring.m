function [t2, spe]=phaseIImonitoring(A,mu0,Sigma0,loadings,THETA)
%
% Function computes the T^2 and SPE charting statistics given in-control
% parameters mu0, Sigma0, loadings, and THETA.
%
% INPUTS:
% A = adjacency matrix
% THE FOLLOWING ARE ESTIMATES OBTAINED FROM PHASE I ANALYSIS
%   mu0 = 16x1 in-control mean vector for triad census
%   Sigma0 = 16x16 in-control variance-covariance matrix for triad census
%   loadings = loadings from Phase I estimation
%   THETA = variance-covariance matrix of the scores 
%
% OUTPUTS:
%   t2 = T-squared charting statistc
%   spe = squared prediction error charting statistic
%
% 
Cc=eye(size(loadings,1))-loadings*loadings';
%
% Compute the linear combinations of triad census proportions
Z1=linearcombs(A,mu0,Sigma0);
%
% compute t's
t=Z1*loadings;
%
% Compute charting statistics, i.e., t2 and spe
t2=t*inv(THETA)*t';
spe=(Z1*Cc)*(Z1*Cc)';