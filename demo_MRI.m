clc;
clear;

addpath('../Common');
addpath('../FRIST_tool');
% load('DataforFiguresandTable/Figure1/I1.mat');
% load('DataforFiguresandTable/Figure1/Q1.mat');
dirdB = 'DataforFiguresandTable';
dirData = {'Table1/Table1_2Drandom5x', 'Table1/Table1_Cartesian7x'};
numData = numel(dirData);

% load('./SaiData/I1.mat');
% load('./SaiData/Q1.mat');
%%%%%%%%%%%%%%%%%%%% MRI parameters %%%%%%%%%%%%%%%%%%%%%%%%%%
% paramsin.nu = 1e6/(512*512);  
paramsin.numiter = 1; 
paramsin.n = 36; 
dim = sqrt(paramsin.n);
% paramsin.N = 512*512;
% paramsin.N = 256 * 256;
paramsin.C = 1e5;

paramsin.lambda0 =0.2;
paramsin.r = 1;
paramsin.W0 = kron(dctmtx(dim),dctmtx(dim));
paramsin.ct = 1;
paramsin.nl = 1;
paramsin.cini = 1;
paramsin.co = 0;
%%%%%%%%%%%%%%%%%%%%% FRIST parameters %%%%%%%%%%%%%%%%%%%%
    FRISTparam.isUnitary = 1;
    FRISTparam.isFlipping = 1;
    FRISTparam.numiterr = 1;
    FRISTparam.mergeK = 40;
    FRISTparam.showStats = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numFRIST = 80;
    constValue = [0.05, 0.055];
    numValue = numel(constValue);

    
    FRIST_result = cell(numData, numValue);
    FRIST_recIm = cell(numData, numValue);
    
    paramsin.stopFlag = 1;
    paramsin.num = numFRIST;
    startValue = 0.012;
    endValue = 0.05;
    numIter = zeros(numValue, numData);
    PSNRpeak = zeros(numValue, numData);
    runTime = zeros(numValue, numData);
    method = {'FRIST'};
    
for idxData = 1 : numData
    for idxValue = 1 : numValue
        load(fullfile(dirdB, dirData{idxData}, 'I1.mat'));
        load(fullfile(dirdB, dirData{idxData}, 'Q1.mat'));

        endValue = constValue(idxValue);
    %     paramsin.s = 0.045*ones(1,paramsin.num);
    %     paramsin.s = [0.012:(0.043/34):0.055  0.055*ones(1, 5)];
        paramsin.s = increasingArray(startValue, endValue, numFRIST - 45);
        paramsin.s = [paramsin.s endValue * ones(1, 45)];
        [FRIST_recIm{idxData, idxValue}, FRIST_result{idxData, idxValue}]= FRIST_MRI(I1,Q1,paramsin, FRISTparam);
        currentPSNR = FRIST_result{idxData, idxValue}.PSNR;
        numIter(idxValue, idxData) = FRIST_result{idxData, idxValue}.numIter;
        PSNRpeak(idxValue, idxData) = max(currentPSNR);
        runTime(idxValue, idxData) = FRIST_result{idxData}.runtime;
        fprintf('FRIST MRI for data = %d has completed!\n', idxData);
        fprintf('FRIST reconstruction PSNR = %.4f.\n', PSNRpeak(idxValue, idxData));

        save('IOP_FRIST_MRI_noisyTuning.mat', 'method', 'PSNRpeak', 'runTime', 'numIter');
    end
end
save('IOP_FRIST_MRI_data_all.mat', 'method', 'FRIST_recIm', 'FRIST_result');

