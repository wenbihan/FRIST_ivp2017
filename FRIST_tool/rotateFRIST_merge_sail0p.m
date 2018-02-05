function [TE, IDX] = rotateFRIST_merge_sail0p(TE, W, permutation, frac, isFlipping, selectIndex)
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
    if ismember(k, selectIndex)
        a1 = TE(L(:,k), :);
        a1 = W * a1;
        a0 = sparseSail0p(a1, frac);
        error(k, :) = sum((a1 - a0).^2);
    else
        error(k, :) = flintmax;
    end
    if isFlipping
        if ismember(k + K, selectIndex)
            a1 = TE(flipInd, :);
            a1 = a1(L(:,k), :);
    %         a1 = TE(L(:,k), :);
    %         a1 = a1(revInd, :);
            a1 = W * a1;
            a0 = sparseSail0p(a1, frac);
            error(k + K, :) = sum((a1 - a0).^2);
        else
            error(k + K, :) = flintmax;
        end
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

end

