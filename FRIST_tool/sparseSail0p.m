function [X] = sparseSail0p(X, frac)
%sparseSail0p Summary of this function goes here
%   Detailed explanation goes here
[s]=sort(abs(X(:)),'descend');
X = X.*(abs(X) >= s(round(frac))); 
end

