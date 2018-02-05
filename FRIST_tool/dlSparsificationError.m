function [RE] = dlSparsificationError(TE, D, CoefMatrix)
%cluster according to the recovered intensityn
RE = (norm(TE - D * CoefMatrix, 'fro')^2) / (norm(TE, 'fro')^2);
end