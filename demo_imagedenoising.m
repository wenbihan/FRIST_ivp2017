clear;
addpath('FRIST_tool');
load('demo_data/noisy_Cameraman_sigma20.mat');

%%%%%%%%%%%%%%%% FRIST image denoising demo %%%%%%%%%%%%%%%%%%%%%
param.n = 64;               % 8 x 8 patches
param.stride = 1;           % fully overlapped
param.sig = sig;           
data.noisy = noisy;
data.oracle = oracle;
clear I7;
param.isKmeansInitialization = false;
param.numBlock = 4;         % 4 rotation directions
[denoised, transform, outputParam]= FRIST_imagedenoising(data, param);
imshow(denoised, []);

