function [flipInd, revInd] = flipIndex_generate(n)
%FLIPINDEX_GENERATE Summary of this function goes here
%   Detailed explanation goes here
indMat = 1 : n;
w = sqrt(n);
indMat = reshape(indMat, [w, w]);
indMat = flip(indMat);
flipInd = indMat(:);
[~, revInd] = sort(flipInd);
end

