function [D, outputParam]= TL_learning(TE, D, param)
%function TL_learning Summary of this function goes here
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs -
%       1. TE : Vectorized data, n * N size 2D matrix.
%       2. D  : Initial transform (example: 2D DCT).
%
%       3. param: Structure that contains the parameters
%                 The various fields are as follows
%       -
%                   - STY: 1 * N row vector, indicates the sparsity levels
%                   - l0: transform Regularizer weight
%                   - iter: number of iterations in learning
%                   - showStats: Set to 1, if the statistics (cost and
%                   sparseError) are output. The run time will be
%                   increased.
%
%
% Outputs -
%       1. D - learned FRIST.
%       2. outputParam: Structure that contains the parameters of the
%       algorithm output for analysis as follows:
%
%                   - time:   run time of the denoising algorithm
%                   (Optional, if showStats = 1)
%                   - sp: sparsifcation error
%
%                   - ct: objective function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
STY = param.STY;
l0 =  param.l0;
isUnitary = param.isUnitary;
iter = param.iter;
showStats = param.showStats;
if ~isUnitary
    l2 = l0 * norm(TE, 'fro')^2;
end
if showStats
    sp = zeros(iter + 1, 1);
    ct = zeros(iter + 1, 1);
    sp(1, 1) = SparsificationError(TE, D, STY);
    if ~isUnitary
        ct(1, 1) = ObjectiveCalc(D, l2, sp(1, 1));
    end
end
tic;
for j =1:iter
    if isUnitary
        D = forwardORTHOTRANSb(D,TE,1,STY);
    else
        D = forwardDL3LearningonlyconjualtcostEXACT2(D, TE, 1, l2,l2, STY);  %param.L is not used in subfunction
    end
    if showStats
        numIter = j + 1;
        sp(numIter, 1) = SparsificationError(TE, D, STY);
        if ~isUnitary
            ct(numIter, 1) = ObjectiveCalc(D, l2, sp(numIter, 1));
        end
    end
end
outputParam.time = toc;
outputParam.ct = ct;
outputParam.sp = sp;
end