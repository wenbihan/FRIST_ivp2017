function [reconstruction] = rotateBack_merge(reconstruction, IDX, permutation)
%ROTATEBACK_MERGE Summary of this function goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs -
%       1. reconstruction : Vectorized data, n * N size 2D matrix.
%       2. IDX:             FRIST index of the rotated data.
%       3. permutation: Structure that contains the permutation
%                 operators
%       -
%                   - L:        permutation vectors n * K
%                   - K:        number of operators
%                   - flipInd:  flipping operator vector
%                   - revInd:   reverse flipping operator vector
%                   - dim:      patch dimension == sqrt(n)
%
% Outputs -
%       1. reconstruction - rotated back version of reconstruction data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = permutation.K;
L = permutation.L;
revInd = permutation.revInd;
isFlipping = permutation.isFlipping;
for k = 1 : K
    [~, rev] = sort(L(:, k));
    reconstruction(:, IDX == k) = reconstruction(rev, IDX == k);
    if isFlipping
        reconstruction(:, IDX == k + K) = reconstruction(rev, IDX == k + K);
        reconstruction(:, IDX == k + K) = reconstruction(revInd, IDX == k + K);
    end
end
end

