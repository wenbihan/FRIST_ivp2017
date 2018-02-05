function [SE] = SparsificationErrorImage(I7, c, n, K, Dc, s, overlap)
%cluster according to the recovered intensity


SE = 0;
[TE0,idx] = my_im2col(I7,[sqrt(n),sqrt(n)],overlap);   %extract patches
mu =  mean(TE0);
TE0 = TE0 - ones(n,1)*mu;
N4 = size(TE0,2);
STY =(ones(1,N4))*(s);      %uniform sparsity constraint

if K==1
    a1 = Dc{1,1}*TE0;
    a0 = sparseSTY(a1, STY);
    SE = (norm(a1-a0,'fro')^2);
else
    a1 = cell(K,1);
a0 = cell(K,1);

for k = 1:K
    a1{k,1} = Dc{k,1}*TE0;
    a0{k,1} = sparseSTY(a1{k,1}, STY);
   SE = SE + (norm(a1{k,1}(:,c{k,1})-a0{k,1}(:,c{k,1}),'fro')^2);
end

end


% if K==1
%     cPse = cell(2,1);
%     for k = 1:2
%         cPse{k,1}=find(IDX==k);
%     a1{k,1} = Dc{1,1}*TE0;
%     a0{k,1} = sparseSTY(a1{k,1}, STY);
%     sizeC = size(cPse{k,1},1);
%     NSE{k,1} = (norm(a1{k,1}(:,cPse{k,1})-a0{k,1}(:,cPse{k,1}),'fro')^2)/(norm(a1{k,1}(:,cPse{k,1}),'fro')^2);
%     wanse = wanse + sizeC*NSE{k,1};
%     end
%     wanse = wanse/N4;
% else