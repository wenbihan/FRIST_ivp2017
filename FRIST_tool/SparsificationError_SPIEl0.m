function [SE] = SparsificationError_SPIEl0(TE, W, threshold)
%cluster according to the recovered intensity
    N = size(TE, 2);
    a1 = W*TE;
    a0 = sparse_SPIEl0(a1, threshold);
    SE = (norm(a1-a0,'fro')^2) / N;
end