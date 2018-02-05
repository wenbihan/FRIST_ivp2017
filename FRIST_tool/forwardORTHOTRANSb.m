function [W,X]= forwardORTHOTRANSb(W,Y,numiter,STY)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs: 1) W : Initial Transform
%        2) Y : Training Matrix with signals as columns
%        3) numiter:  Number of iterations of alternating minimization
%        4) STY: Vector containing maximum allowed sparsity levels for each training signal.

%Outputs:  1) W: Learnt Transform
%          2) X: Learnt Sparse Code


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial steps
[K,n]=size(W);

ix=find(STY>0); q=Y(:,ix); STY=STY(:,ix); N=size(q,2);
ez=K*(0:(N-1));STY=STY + ez;

%Algorithm iterations in a FOR Loop
for i=1:numiter

    %Sparse Coding Step
    X1=W*q;
    [s]=sort(abs(X1),'descend');
    X = X1.*(bsxfun(@ge,abs(X1),s(STY)));
    
    %Transform Update Step
    [U,~,V]=svd(q*X');
    W=(V(:,1:n))*U';

end