function [X] = sparse_SPIEl0(X, frac)
%sparse_SPIEl0 Summary of this function goes here
%   Detailed explanation goes here
X = X.*(abs(X) >= frac);
end

