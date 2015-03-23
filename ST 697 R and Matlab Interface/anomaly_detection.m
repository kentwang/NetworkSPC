function [loadings,THETA,UCLt,UCLs]=anomaly_detection(X,alpha)
%
% This function impements the methods outlines in Alberto Ferrer's 2014
% Quality Engineering paper on monitoring Multivariate processes.
%
% Ferrer, A. (2014). "Latent Structures-Based Mutivariate Statistical
% Process Control: A Paradigm Shift", Quality Engineering 26:1 72-91.
%
% Inputs: 
% (1) X = m x 9 matrix of in-control vectors of linear combinations of triad census 
% (2) alpha = type I error probability

% Standardize input observations
for k=1:size(X,2)
    Z(:,k)=(X(:,k)-mean(X(:,k)))/sqrt(var(X(:,k)));
end
%
% Compute covariance matrix of Z 
R=cov(Z);
% Note: since Z is standardized, R is correlation matrix
%
% Compute eigenvalues/eigenvectors of correlation matrix R
[V,D]=eig(R);
d=diag(D);
d=sort(d,'descend');

% Retain components using eigen-one criteria
j=0;
for i=1:length(d)
    if d(i)>=1
        j=j+1;
        loadings(:,j)=V(:,length(d)-i+1);
    end
end

U=Z*V;

% Calculate scores
scores=Z*loadings;

% Calculate variance-covariance matrix of scores
THETA=cov(scores);

% Calculate residuals
E=Z-scores*loadings';

% Generate T2 and SPE statistics
for i=1:size(X,1)
    % T-squared statistic
    T2(i)=scores(i,:)*cov(scores)^(-1)*scores(i,:)';
    % squared prediction error statistic
    SPE(i)=E(i,:)*E(i,:)';
end

% Phase I control limit for T2 control chart
m=size(X,1);
CL1=((m-1)^2/m)*betainv(1-alpha,size(scores,2)/2,(m-size(scores,2)-1)/2);
%
% Phase II control limit for T2 control chart
UCLt=(size(scores,2)*(m^2-1)/(m*(m-size(scores,2))))*finv(1-alpha,size(scores,2),m-size(scores,2));

% Control limit for SPE control chart
phi1=0; phi2=0; phi3=0;
for k=j+1:length(d)
    phi1=phi1+d(k)^1;
    phi2=phi2+d(k)^2;
    phi3=phi3+d(k)^3;
end
h0=1-2*phi1*phi3/(3*phi2^2);
CL2=phi1*(norminv(1-alpha,0,1)*sqrt(2*phi2*h0^2)/phi1 + 1 + phi2*h0*(h0-1)/phi1^2)^(1/h0);
UCLs=CL2;
%
%
%Plot Phase I control chart
for i=1:m
    Limit(i)=CL1;
    Limit2(i)=CL2;
end
%
i=1:m;
subplot(2,1,1); plot(i,T2,'ko--',i,Limit,'r'); title('T^2 Chart on Scores'); ylabel('T^2');
subplot(2,1,2); plot(i,SPE,'ko--',i,Limit2,'r'); title('SPE Chart on Residual Distance'); xlabel('Time'); ylabel('SPE');
% 
% Probability of signal
% p1=sum(T2>CL1)/m;
% p2=sum(SPE>CL2)/m;
% 
% Find indices of out of control points
% I1=find(T2>CL1);
% I2=find(SPE>CL2);

% Xn and Xa are the observations projected onto the normal and anomolous subspaces
% Note: X = Xn + Xa
%
% Normal subspace projection
%Xn=X*loadings*loadings';
%
% Anomolous subspace projection (i.e., residuals)
%Xa=X*(eye(size(loadings,1))-loadings*loadings');
