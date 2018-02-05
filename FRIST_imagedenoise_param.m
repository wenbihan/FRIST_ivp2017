function [param] = FRIST_imagedenoise_param(param)
%Function for tuning paramters for FRIST denoising
% Input:
%       - n         :  patch size
%       - stride    :  patch stride in extraction
%       - sig       :  noise standard deviation
%       - numBlock  :  #directions & rotations
%       - isKmeansInitialization    :   use Kmeans initialization?
%Note that all input parameters need to be set prior to simulation. 
%This tuning function is just an example settings which we provide, for 
%generating the results in the "FRIST paper". However, the user is
%advised to carefully modify this function, thus choose optimal values  
%for the parameters depending on the specific data or task at hand.

% options
param.isFlipping = 1;
param.isUnitary = 0;
param.isMultipass = true;
n   =   param.n;

param.mergeK    =   param.numBlock;
sig = param.sig;
isMultipass = param.isMultipass;
%% constant     -   C1
if sqrt(n) <= 8
    C1 = 1.08;
elseif sqrt(n) <= 9
    C1 = 1.07;
elseif sqrt(n) <= 10
    C1 = 1.05;
elseif sqrt(n) <= 11
    C1 = 1.04;
elseif sqrt(n) <= 12
    C1 = 1.03;
elseif sqrt(n) <= 15
    C1 = 1.015;
elseif sqrt(n) <=20
    C1 = 1.01;
else
    C1 = 1.0;
end
%% multi-pass denoising signal list     -   sig2
if isMultipass
    if sig <= 8
        sig2 = sig;
    elseif sig <= 12
        sig2 = [sig*0.9, sig/3];
    elseif sig <= 17
        sig2 = [sig*0.9, sig*0.2, sig*0.1, sig*0.05];
    elseif sig <= 50
        sig2 = [sig*0.9, sig*0.2, sig*0.1];
    else
        sig2 = [sig*0.9, sig*0.2, sig*0.04, sig*0.03];
    end
else
    sig2 = sig;
end
%% coefficient  -   la
la = 0.01/sig;
%% initial sparsity     -   T0
T0 = round((0.1)*n);
%% initial OCTOBOS   -   transform
% Transform Initialization: DCT
D = kron(dctmtx(sqrt(n)),dctmtx(sqrt(n)));  
if isfield(param, 'isOTL') && param.isOTL == 1
    D = [D; eye(n)];
    D = D(1:100, :);
end
        
%% number of iteration    -   iter
iter = 40;
param.C1 = C1;
param.sig2 = sig2;
param.la = la;
param.T0 = T0;
param.transform = D;
param.iter = iter;
param.l0 = 0.031;
% param.iterMultipass = 3;
param.iterMultipass = 20;
% param.roundLearning = 12;
param.roundLearning = 10;

end

