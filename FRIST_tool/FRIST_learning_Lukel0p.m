function [YH, D, outputParam]= FRIST_learning_Lukel0p(TE, D, param)
%function FRIST_learning_merge Summary of this function goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs -
%       1. TE : Vectorized data, n * N size 2D matrix.
%       2. D  : Initial transform (example: 2D DCT).
%
%       3. param: Structure that contains the parameters
%                 The various fields are as follows
%       -
%                   - isFlipping: set to 1, if flipping is considered
%                   - numiterr:  number of intra iterations (between sparse
%                   coding and transform update)
%                   - numIt:      number of outer iterations
%                   - showStats: Set to 1, if the statistics (cost and
%                   sparseError) are output. The run time will be
%                   increased.
%                   - mergeK:    number of desired number of clusters
%                   - thresh:      threshold for l0 penalty
%
%
% Outputs -
%       1. YH- permuted version of TE 
%       2. D - learned FRIST.
%       3. outputParam: Structure that contains the parameters of the
%       algorithm output for analysis as follows:
%
%                   - IDX :   clustering result indices
%                   - permutation: permutation parameters.
%                   (Optional, if showStats = 1)
%                   - sp: sparsifcation error
%                   - numCl :  numbers of cluster element changes
%                   - selectIndex: selected operator index
%                   - permutation: generated permutation struct.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
threshold = param.thresh;
% l0 =  param.l0;
% if isfield(param, 'isUnitary')
%     isUnitary = param.isUnitary;
% end
iter = param.numIt;
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
% if isfield(param, 'isUnitary') && ~isUnitary
%     l2 = l0 * norm(TE, 'fro')^2;
% end
if isfield(param, 'isFlipping') && isFlipping
    numCl = zeros(2 * permutation.K, iter);
    if isfield(param, 'giveSelectIndex') && param.giveSelectIndex
        selectIndex = param.selectIndex;
        reduceK = numel(selectIndex);
    else
        selectIndex = 1 : (2 * K);
        reduceK = 2 * K;
    end
else
    numCl = zeros(permutation.K, iter);
    if isfield(param, 'giveSelectIndex') && param.giveSelectIndex
        selectIndex = param.selectIndex;
        reduceK = numel(selectIndex);
    else
        selectIndex = 1 : K;
        reduceK = K;
    end

end
if showStats
    sp = zeros(iter + 1, 1);
    sp(1, 1) = SparsificationError_SPIEl0(TE, D, threshold);
%     if ~isUnitary
%         ct = zeros(iter + 1, 1);
%         ct(1, 1) = ObjectiveCalc(D, l2, sp(1, 1));
%     end
end
% tic;
for j =1:iter
    %%%%%%%%%%% clustering %%%%%%%%%%
    [YH, IDX] = rotateFRIST_merge_SPIEl0(TE, D, permutation, threshold, isFlipping, selectIndex);
    %%%%%%%%%%% transform learning over FR data %%%%%%%%%
    numiter = numiterr;
%     if isfield(param, 'isUnitary') && isUnitary
%         D = forwardORTHOTRANSb_SPIEl0(D,YH,numiter,threshold);
%     else
%         D = forwardDL3LearningonlyconjualtcostEXACT2_SPIEl0(D, YH, numiter, l2,l2, threshold);  %param.L is not used in subfunction
%     end
    D = forwardORTHOTRANSb_SPIEl0(D,YH,numiter,threshold);
    if showStats
        numIter = j + 1;
        sp(numIter, 1) = SparsificationError_SPIEl0(YH, D, threshold);
%         if ~isUnitary
%             ct(numIter, 1) = ObjectiveCalc(D, l2, sp(numIter, 1));
%         end
    end
    for k = 1 : K        
        numCl(k, j) = numel(find(IDX == k));
        if isFlipping
            numCl(k + K, j) = numel(find(IDX == k + K));
        end
    end
    %%%%%%%%% select the largest K cluster %%%%%%%%%%%%%
    while reduceK > mergeK
        reduceK = max(round(reduceK / 2), mergeK);
        [~, sortweight] = sort(numCl(:, j),'descend');
        selectIndex = sortweight(1 : reduceK);
    end
%     if j == 3
%         [~, sortweight] = sort(numCl(:, j),'descend');
%         selectIndex = sortweight(1 : mergeK);
%     end
end
% outputParam.time = toc;
outputParam.IDX = IDX;
if showStats
    outputParam.sp = sp;
%     if ~isUnitary
%         outputParam.ct = ct;
%     end
    outputParam.numCl = numCl;
end 
outputParam.selectIndex = selectIndex;
outputParam.permutation = permutation;
outputParam.reduceK = reduceK;
end