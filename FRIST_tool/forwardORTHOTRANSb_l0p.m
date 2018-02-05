function [W,X]= forwardORTHOTRANSb_l0p(W,Y,numiter,l5)

%This is an implementation of the orthonormal or unitary transform learning algorithm with closed-form solutions for the sparse coding and transform update steps that was presented in the following papers:
%1) S. Ravishankar and Y. Bresler, "Closed-form solutions within sparsifying transform learning," in Proc. IEEE International Conference on Acoustics, Speech and Signal Processing (ICASSP), 2013, pp. 5378-5382.
%2) S. Ravishankar and Y. Bresler, “Closed-Form Optimal Updates In Transform Learning,?in Signal Processing with Adaptive Sparse Structured Representations (SPARS) workshop, July 2013.
%3) S. Ravishankar and Y. Bresler, “\ell_0 Sparsifying transform learning with efficient optimal updates and convergence guarantees,?IEEE Transactions on Signal Processing, vol. 63, no. 9, pp. 2389-2404, May 2015.

%We employ alternating minimization here to solve a transform learning problem that involves a constraint on the adaptive transform domain sparsity of each training signal, 
%and a constraint that the transform is unitary or orthonormal.
%The algorithm iterates over a sparse coding step and a transform update step, both of which involve efficient update procedures.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Inputs: 1) W : Initial Transform
%        2) Y : Training Matrix with signals as columns
%        3) numiter:  Number of iterations of alternating minimization
%        4) STY: Vector containing maximum allowed sparsity levels for each training signal.

%Outputs:  1) W: Learnt Transform
%          2) X: Learnt Sparse Code


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initial steps
[~,n]=size(W);
% N = size(Y, 2);

% ix=find(STY>0); q=Y(:,ix); STY=STY(:,ix); N=size(q,2);
% ez=K*(0:(N-1));STY=STY + ez;

%Algorithm iterations in a FOR Loop
for i=1:numiter

    %Sparse Coding Step
    X1 = W * Y;
%     [s] = sort(abs(X1),'descend');
    st= (l5)*ones(1,size(X1,2));    
    X = X1.*(bsxfun(@ge,abs(X1),st));
    
    %Transform Update Step
    [U,~,V]=svd(Y*X');
    W=(V(:,1:n))*U';
end