clear;
addpath('FRIST_tool');
load('demo_data/I1');
load('demo_data/Q1');

%%%%%%%%%%%%%%%% FRIST MRI demo %%%%%%%%%%%%%%%%%%%%%
paramsin.numiter = 1; 
paramsin.n = 36; 
dim = sqrt(paramsin.n);
paramsin.C = 1e5;
paramsin.lambda0 =0.2;
paramsin.r = 1;
paramsin.W0 = kron(dctmtx(dim),dctmtx(dim));
paramsin.ct = 1;
paramsin.nl = 1;
paramsin.cini = 1;
paramsin.co = 0;
paramsin.num = 80;
paramsin.stopFlag = 1;
FRISTparam.isUnitary = 1;
FRISTparam.isFlipping = 1;
FRISTparam.numiterr = 1;
FRISTparam.mergeK = 40;
FRISTparam.showStats = 1;
startValue  =   0.012;
endValue    =   0.05;
numFRIST    =   80;
paramsin.s  =   increasingArray(startValue, endValue, numFRIST - 45);
paramsin.s  =   [paramsin.s endValue * ones(1, 45)];
[Xr, outputParam] = ...
    FRIST_MRI(I1, Q1, paramsin, FRISTparam);
% imshow(Xr, []);

