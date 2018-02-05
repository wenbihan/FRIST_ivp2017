function [Xr, Dlist, outputParam]= FRIST_imagedenoising(data, param)
%Function for denoising the gray-scale image using FRIST-based denoising
%algorithm.
%
%Note that all input parameters need to be set prior to simulation. We
%provide some example settings using function FRIST_imagedenoisinge_param.
%However, the user is advised to carefully choose optimal values for the
%parameters depending on the specific data or task at hand.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Inputs -
%       1. data : Image data. The fields are as follows -
%                   - noisy: a*b size gray-scale matrix for denoising
%                   - oracle (optional): a*b size gray-scale matrix as
%                   ground-true, which is used to calculate PSNR
%
%       2. param: Structure that contains the parameters of the
%       OCTOBOS_imagedenoising algorithm. The various fields are as follows
%       -
%                   - isFlipping: whether including flipping
%                   - sig: Standard deviation of the additive Gaussian
%                   noise (Example: 20)
%                   - n: Patch size as (Example: 64)
%                   - stride: stride of overlapping patches
%                   - isUnitary: set to 1, if unitary TL is used
%                   - isMultipass: set to 1, if multipass is applied
%
% Outputs -
%       1. Xr - Image reconstructed with FRIST_imagedenoising algorithm.
%       2. Dlist - learned FRIST.
%       2. outputParam: Structure that contains the parameters of the
%       algorithm output for analysis as follows
%       -
%                   - psnrOut: PSNR of Xr, if the oracle is provided
%                   - time:   run time of the denoising algorithm

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%  Initialization  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param = FRIST_imagedenoise_param(param);
mergeK = param.mergeK;
% n = param.n;                                % patch size / dimensionality
isUnitary = param.isUnitary;
C1 = param.C1;                              % thresholding coefficient
sig2 = param.sig2;                          % multi-pass noise level estimates
la = param.la;                              % fidelity term coefficient
T0 = param.T0;                              % initial sparsity level
Xr = data.noisy;                            % noisy image
transform = param.transform;                % initial transform
iter = param.iter;                          % number of iterations in first pass denoising
iterMultipass = param.iterMultipass;        % number of iterations in multipass denoising
l0 = param.l0;                              % regularizer coefficient
roundLearning = param.roundLearning;        % number of rounds for OCTOBOS learning
stride = param.stride;                      % stride of overlapping patches
isFlipping = param.isFlipping;
% clear param;
stp = 1;                                    % sparsity increase stepsize
SP = 1:stp:round(9*T0);                     % maximum sparsity level allowed in algorithm is 9*T0 here
len = length(sig2);
[nK, n] = size(transform);
Dlist = zeros(nK, n, len);
permutation = generatePermutation(sqrt(n), isFlipping);
L = permutation.L;
K = permutation.K;
revInd = permutation.revInd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  MAIN CODE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
for pr = 1 : length(sig2)
    if isFlipping
        numCl = zeros(2 * permutation.K, iter);
        selectIndex = 1 : (2 * K);
        reduceK = 2 * K;
    else
        numCl = zeros(permutation.K, iter);
        selectIndex = 1 : K;
        reduceK = K;
    end
    if(pr > 1)
        iter = iterMultipass;
    end
    sig = sig2(pr);
    threshold = C1 * sig * (sqrt(n));                   %threshold in variable sparsity update    
    [TE, idx] = my_im2col(Xr,[sqrt(n),sqrt(n)],stride);
    NTE = size(TE,2);
    mu = mean(TE);
    TE = TE - ones(n,1)*mu;                           %mean subtraction
    STY = ones(1, NTE).*T0;                            % Initial Sparsity Vector
    D = transform;
    if ~isUnitary
        l2 = l0 * norm(TE, 'fro')^2;
    end
    for j =1 : iter
        [YH, IDX] = rotateFRIST_merge(TE, D, permutation, STY, isFlipping, selectIndex);
        if j < iter/3
            currentRoundLearning = roundLearning*2;    %roundLearning is maximum number of learning iterations
        else
            currentRoundLearning = roundLearning;
        end
        if isUnitary
            % unitary transform learning
            D = forwardORTHOTRANSb(D, YH,currentRoundLearning, STY);
        else
            % well-conditioned transform learning
            [D, ~] = transformLearning(YH, D, l2, l2, STY, currentRoundLearning);
        end
        [STY, reconstruction] = sparsityUpdate(YH, D, la, threshold, SP);        % Sparsity Update for YH
        %%%%%%%%% select the largest K cluster %%%%%%%%%%%%%
        while reduceK > mergeK
            for k = 1 : K
                numCl(k, j) = numel(find(IDX == k));
                if isFlipping
                    numCl(k + K, j) = numel(find(IDX == k + K));
                end
            end
            reduceK = max(round(reduceK / 2), mergeK);
            [~, sortweight] = sort(numCl(:, j),'descend');
            selectIndex = sortweight(1 : reduceK);
        end
    end
%  Patches Recovery
    for k = 1 : K
        [~, rev] = sort(L(:, k));
        reconstruction(:, IDX == k) = reconstruction(rev, IDX == k);
        if isFlipping
            reconstruction(:, IDX == k + K) = reconstruction(rev, IDX == k + K);
            reconstruction(:, IDX == k + K) = reconstruction(revInd, IDX == k + K);
        end
    end
    reconstruction = reconstruction + ones(n,1)*mu;
    Xr = my_col2image(Xr, n, reconstruction, idx, 0);
    Xr(Xr>255) = 255;
    Xr(Xr<0) = 0;
    Dlist(:, :, pr) = D;
%     display(pr);
end
outputParam.time = toc;
outputParam.psnrOut = PSNR(Xr - data.oracle);
end