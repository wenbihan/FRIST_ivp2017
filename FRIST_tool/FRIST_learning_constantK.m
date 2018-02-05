function [YH, D, outputParam]= FRIST_learning_constantK(TE, D, param)
%function FRIST_learning Summary of this function goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs -
%       1. TE : Vectorized data, n * N size 2D matrix.
%       2. D  : Initial transform (example: 2D DCT).
%
%       3. param: Structure that contains the parameters
%                 The various fields are as follows
%       -
%                   - STY:       1 * N row vector, indicates the sparsity levels
%                   - l0:        transform Regularizer weight
%                   - isUnitary: set to 1, if applying unitary transform
%                   - isFlipping: set to 1, if flipping is considered
%                   - numiterr:  number of intra iterations (between sparse
%                   coding and transform update)
%                   - iter:      number of outer iterations
%                   - showStats: Set to 1, if the statistics (cost and
%                   sparseError) are output. The run time will be
%                   increased.
%
%
% Outputs -
%       1. YH- permuted version of TE 
%       2. D - learned FRIST.
%       3. outputParam: Structure that contains the parameters of the
%       algorithm output for analysis as follows:
%
%                   - time:   run time of the denoising algorithm
%                   - IDX :   clustering result indices
%                   - permutation: permutation parameters.
%                   (Optional, if showStats = 1)
%                   - sp: sparsifcation error
%
%                   - ct: objective function
%                   - numCl :  numbers of cluster element changes

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STY = param.STY;
l0 =  param.l0;
isUnitary = param.isUnitary;
iter = param.iter;
showStats = param.showStats;
isFlipping = param.isFlipping;
numiterr = param.numiterr;
dim = sqrt(size(TE, 1));
permutation = generatePermutation(dim, isFlipping);
%%%%%%%%% manually select the operators from all available ones %%%%%%
selectIndex = param.selectIndex;
permutation.K = numel(selectIndex);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
permutation.L = permutation.L(:, selectIndex);
K = permutation.K;
% L = permutation.L;
% if isFlipping
%     flipInd = permutation.flipInd;
%     revInd = permutation.revInd;
% end
if ~isUnitary
    l2 = l0 * norm(TE, 'fro')^2;
end
if showStats
    Dchange = zeros(iter, 1);
    sp = zeros(iter + 1, 1);
    sp(1, 1) = SparsificationError(TE, D, STY);
    if isFlipping
        numCl = zeros(2 * permutation.K, iter);
    else
        numCl = zeros(permutation.K, iter);
    end
    if ~isUnitary
        ct = zeros(iter + 1, 1);
        ct(1, 1) = ObjectiveCalc(D, l2, sp(1, 1));
    end
end
tic;
for j =1:iter
%%%%%%%%%%% clustering %%%%%%%%%%
[YH, IDX] = rotateFRIST(TE, D, permutation, STY, isFlipping);
%%%%%%%%%%% transform learning over FR data %%%%%%%%%
    numiter = numiterr;
%     if(j<(iter/3))
%         numiter=numiterr*2;    %numiterr is maximum number of learning iterations
%     else
%         numiter=numiterr;
%     end
    D1 = D;
    if isUnitary
        D = forwardORTHOTRANSb(D,YH,numiter,STY);
    else
        D = forwardDL3LearningonlyconjualtcostEXACT2(D, YH, numiter, l2,l2, STY);  %param.L is not used in subfunction
    end
    if showStats
        numIter = j + 1;
        sp(numIter, 1) = calSparsificationError(YH, D, STY);
        Dchange(j, 1) = norm(D - D1, 'fro');
        if ~isUnitary
            ct(numIter, 1) = ObjectiveCalc(D, l2, sp(numIter, 1));
        end
    end
    for k = 1 : K
        numCl(k, j) = numel(find(IDX == k));
        if isFlipping
            numCl(k + K, j) = numel(find(IDX == k + K));
        end
    end
end
outputParam.time = toc;
outputParam.IDX = IDX;
outputParam.Dchange = Dchange;
if showStats
    outputParam.sp = sp;
    if ~isUnitary
        outputParam.ct = ct;
    end
    outputParam.numCl = numCl;
end 
outputParam.permutation = permutation;
outputParam.regularizer = l2 * (norm(D, 'fro')^2 - log(abs(det(D))));
end