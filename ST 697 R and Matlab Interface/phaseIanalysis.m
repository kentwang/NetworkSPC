function [mu, Sigma, loadings, theta, uclT2, uclSPE]=phaseIanalysis(adj,alpha)
%
% This function performs a Phase I analysis of a series of directed graphs
% using the method suggested by Alberto Ferrer.
%
% Inputs: 
%   adj = a series of adjacency matrices stored in a 1xm cell array
%   alpha = type I error probability
%
% Outputs: 
%   Needed to compute the charting statistics in Phase II:
%       mu = 16x1 in-control mean vector of the triad census proportions
%       Sigma = 16x16 in-control covariance matrix of triad census proportions
%       loadings = 9xk loadings matrix from the PCA
%       theta = kxk Covariance matrix of the scores from the PCA 
%       (Note: k denotes the number of principle components retained)
%
%   Control limits for use in subsequent Phase II:
%       uclT2 = Phase II control limit for T-squared chart
%       uclSPE = Phase II control limit for SPE chart
%
%
%
% Compute triad census for each of the m in-control historical networks
for i=1:length(adj)
    T(i,1:16)=triadcensus(adj{i});
end
%
% Compute mean and variance-covariance matrix of the triad census
mu=mean(T)'; Sigma=cov(T);
%
% Compute linear combinations of the triad census
for i=1:length(adj)
    Z(i,1:9)=linearcombs(adj{i},mu,Sigma);
end
%
% Perform phase I analysis using methods proposed in Ferrer's 2014 QE paper
[loadings,theta,uclT2,uclSPE]=anomaly_detection(Z,alpha);
