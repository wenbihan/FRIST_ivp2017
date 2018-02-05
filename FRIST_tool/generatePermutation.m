function [permutation] = generatePermutation(dim, isFlipping)
%GENERATEPERMUTATION Summary of this function goes here
%   Detailed explanation goes here
%%%% generate rotational permutations %%%%%
n = dim^2;
L = compute_geometric_lut_limitedAngle(dim);
K = size(L, 2);
permutation.L = L;
%%%% generate flipping permutations %%%%%
if isFlipping
    [flipInd, revInd] = flipIndex_generate(n);
    permutation.flipInd = flipInd;
    permutation.revInd = revInd;
end
permutation.K = K;
permutation.dim = dim;
permutation.isFlipping = isFlipping;
end

