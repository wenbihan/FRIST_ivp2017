function [SE] = SparsificationError_sail0p(TE, W, frac)
%cluster according to the recovered intensity
    N = size(TE, 2);
    a1 = W*TE;
    a0 = sparseSail0p(a1, frac);
    SE = (norm(a1-a0,'fro')^2) / N;
end