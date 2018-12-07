%% -- DISM with brain data
clear all;close all; clc;
load data/datSDM; 

data = dataAll{3};
% data.W = data.Wnorm;

data.X = Xnorm(data.X,1);
opt.lambda1s = [logspace(2.5, 0, 9) 0]/1000;
opt.lambda2s = [0 0.1 0.5 1 2];

opt.lambda1s = [0.15];
opt.lambda2s = [0.1];


opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
% supLE = supLE(data,opt);
tStart = tic; m = DISM(data, opt); toc(tStart) % m.accuCV

mSVM = accuracySVM(data);

disp(strcat('avg svm accuracy:',num2str(mean(mSVM.accuCV)), '(' ,num2str(std(mSVM.accuCV)), ')'));

visual(data, m);



% Experiment 1: fix lambda1 = 0.07 and k = 50, change lambda2, increase the
% number of selected features and visualize selected subnetworks



% Experiment 2: fix k = 50, change lambda1 and lambda2 at the same time,
% fix the number of selected features and visualize selected subnetworks