function T=triadcensus(A)
%
% This function computes the triad census as given in Holland and Leinart's
% original paper.  Both directed and undirected networks are allowed.
% 
% Inputs: adjacency matrix A
% Outputs: Triad census proportions T
%   
n=size(A,1);
%
E=zeros(n,n);
M=zeros(n,n);
for i=1:n
    for j=1:n
        if A(i,j)==1||A(j,i)==1
            E(i,j)=1; E(j,i)=1;
        end
        if A(i,j)==1 && A(j,i)==1
            M(i,j)=1; M(j,i)=1;
        end
    end
end
%
% Define the following for use in counting configurations
C=A-M; Ec=ones(n,n)-eye(n)-E;
%
% Counts of 16 configurations in the triad census
T003=trace(Ec^3)/6;
T012=sum(sum(Ec^2.*(C+C')))/2;
T102=sum(sum(Ec^2.*M))/2;
T021D=sum(sum(C'*C.*Ec))/2;
T021U=sum(sum(C*C'.*Ec))/2;
T021C=sum(sum(C^2.*Ec));
T111D=sum(sum(A*A'.*Ec-M^2.*Ec-C*C'.*Ec))/2;
T111U=sum(sum(A'*A.*Ec-M^2.*Ec-C'*C.*Ec))/2;
T030T=sum(sum(C^2.*C));
T030C=trace(C^3)/3;
T201=sum(sum(M^2.*Ec))/2;
T120D=sum(sum(C'*C.*M))/2;
T120U=sum(sum(C*C'.*M))/2;
T120C=sum(sum(C^2.*M));
T210=sum(sum(M^2.*(C+C')))/2;
T300=trace(M^3)/6;
%
% Triad census counts
T=[T003 T012 T102 T021D T021U T021C T111D T111U T030T T030C T201 T120D T120U T120C T210 T300];
%
% Triad census proportions
T=T/nchoosek(size(A,1),3);
