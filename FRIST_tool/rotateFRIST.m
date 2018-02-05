function [TE, IDX] = rotateFRIST(TE, W, permutation, STY, isFlipping)
%SPARSESTY Summary of this function goes here
%   Detailed explanation goes here
K = permutation.K;
NTE = size(TE, 2);
L = permutation.L;
if isFlipping
    flipInd = permutation.flipInd;
    revInd = permutation.revInd;
    error = zeros(K * 2, NTE);
else
    error = zeros(K, NTE);
end
for k = 1 : K
    a1 = TE(L(:,k), :);
    a1 = W * a1;
    a0 = sparseSTY(a1, STY);
    error(k, :) = sum((a1 - a0).^2);
    if isFlipping
        a1 = TE(flipInd, :);
        a1 = a1(L(:,k), :);
%         a1 = TE(L(:,k), :);
%         a1 = a1(revInd, :);
        a1 = W * a1;
        a0 = sparseSTY(a1, STY);
        error(k + K, :) = sum((a1 - a0).^2);
    end
end
[~, IDX]=min(error, [], 1);
IDX = IDX';
for k = 1 : K
    TE(:, IDX == k) = TE(L(:, k), IDX == k);
    if isFlipping
        TE(:, IDX == k + K) = TE(flipInd, IDX == k + K);
        TE(:, IDX == k + K) = TE(L(:, k), IDX == k + K);
    end
end

% [K, N] = size(X1);
% ix=find(STY>0); 
% a0=X1(:,ix); 
% STY=STY(:,ix); 
% N=size(a0,2);
% ez=K*(0:(N-1));
% STY = STY + ez;
% [s]=sort(abs(a0),'descend');
% a1 = a0.*(bsxfun(@ge,abs(a0),s(STY)));
% X = zeros(size(X1));
% X(:,ix) = a1;
end

