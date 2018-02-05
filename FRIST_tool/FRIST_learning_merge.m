function [YH, D, outputParam]= FRIST_learning_merge(TE, D, param)
%function FRIST_learning_merge Summary of this function goes here
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
%                   - mergeK:    number of desired number of clusters
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
%                   - selectIndex: selected operator index

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STY = param.STY;
l0 =  param.l0;
if isfield(param, 'isUnitary')
    isUnitary = param.isUnitary;
end
if isfield(param, 'isOTL')
    isOTL = param.isOTL;
end
iter = param.iter;
if isfield(param, 'showStats')
    showStats = param.showStats;
end
if isfield(param, 'isFlipping')
    isFlipping = param.isFlipping;
end
numiterr = param.numiterr;
dim = sqrt(size(TE, 1));
permutation = generatePermutation(dim, isFlipping);
K = permutation.K;
mergeK = param.mergeK;
% L = permutation.L;
% if isFlipping
%     flipInd = permutation.flipInd;
%     revInd = permutation.revInd;
% end
if isfield(param, 'isOTL') && ~isUnitary
    l2 = l0 * norm(TE, 'fro')^2;
end
if isfield(param, 'isFlipping') && isFlipping
    numCl = zeros(2 * permutation.K, iter);
    selectIndex = 1 : (2 * K);
    reduceK = 2 * K;
else
    numCl = zeros(permutation.K, iter);
    selectIndex = 1 : K;
    reduceK = K;
end
if showStats
    sp = zeros(iter + 1, 1);
    sp(1, 1) = SparsificationError(TE, D, STY);
    if ~isUnitary
        ct = zeros(iter + 1, 1);
        ct(1, 1) = ObjectiveCalc(D, l2, sp(1, 1));
    end
end
tic;
counter = 0;
for j =1:iter
%%%%%%%%%%% clustering %%%%%%%%%%
[YH, IDX] = rotateFRIST_merge(TE, D, permutation, STY, isFlipping, selectIndex);
%%%%%%%%%%% transform learning over FR data %%%%%%%%%
    numiter = numiterr;
%     if(j<(iter/3))
%         numiter=numiterr*2;    %numiterr is maximum number of learning iterations
%     else
%         numiter=numiterr;
%     end
    if isfield(param, 'isUnitary') && isUnitary
        D = forwardORTHOTRANSb(D,YH,numiter,STY);
    else
        if isOTL
            mu = param.mu;
            numgrad = param.numgrad;
            p = param.p;
            D = TLTALL(D, YH, numiter, l2,l2,l2, p,mu,numgrad,STY);
        else
            D = forwardDL3LearningonlyconjualtcostEXACT2(D, YH, numiter, l2,l2, STY);  %param.L is not used in subfunction
        end
    end
    if showStats
        numIter = j + 1;
        sp(numIter, 1) = SparsificationError(YH, D, STY);
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
    %%%%%%%%% select the largest K cluster %%%%%%%%%%%%%
    while reduceK > mergeK
        counter = counter + 1;
        reduceK = max(round(reduceK / 2), mergeK);
        [~, sortweight] = sort(numCl(:, j),'descend');
        selectIndex = sortweight(1 : reduceK);
    end
%     if j == 3
%         [~, sortweight] = sort(numCl(:, j),'descend');
%         selectIndex = sortweight(1 : mergeK);
%     end
end
outputParam.time = toc;
outputParam.IDX = IDX;
if showStats
    outputParam.sp = sp;
    if ~isUnitary
        outputParam.ct = ct;
    end
    outputParam.numCl = numCl;
end 
outputParam.permutation = permutation;
outputParam.selectIndex = selectIndex;
end