% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 1, CMU straight-poses, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% DISM
clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'Exp-1-CMU-all-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = true;
opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
% opt.lambda1s = 0.050:-0.001:0.045;
% opt.lambda1s = 0.025:-0.002:0.005;
opt.lambda1s =  [0.0820, 0.0815, 0.0810, 0.0805];
opt.lambda2s = 350; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);



% DISMw/o
clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'Exp-1-CMU-straight-DISMwo';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = true;
opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s =  0.02;
opt.lambda2s = 0; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);


% PCA SVM
[U, S, V]=SVD_CutOff(data.X,.95);
data.X = U'*data.X;
mSVM = accuracySVM(data); 
disp(sprintf('acc: %6.4f, std: %6.4f', mean(mSVM.accuCV), std(mSVM.accuCV)));


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 7, CMU straight-poses, AUC values
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'Exp-7-CMU-straight-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

% DISM w/o
disp('DISM w/o');
opt.bLinear = 0; opt.alpha = .2; opt.output = false;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda2s = 0; opt.k = 50;
% 25 features
opt.lambda1s = 0.065; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 50 features
opt.lambda1s = 0.04; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 100 features
opt.lambda1s = 0.01; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 200 features
opt.lambda1s = 0.0004; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));



% DISM
disp('DISM');
opt.bLinear = 0; opt.alpha = .2; opt.output = false;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda2s = 350; opt.k = 50;
% 25 features
opt.lambda1s = 0.15; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 50 features
opt.lambda1s = 0.105; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 100 features
opt.lambda1s = 0.075; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 200 features
opt.lambda1s = 0.055; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 7, CMU all-poses, AUC values
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{2};

description = 'Exp-7-CMU-all-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

% DISM w/o
disp('DISM w/o');
opt.bLinear = 0; opt.alpha = .2; opt.output = false;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda2s = 0; opt.k = 50;
% 25 features
opt.lambda1s = 0.09; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 50 features
opt.lambda1s = 0.065; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 100 features
opt.lambda1s = 0.04; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 200 features
opt.lambda1s = 0.0225; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));



% DISM
disp('DISM');
opt.bLinear = 0; opt.alpha = .2; opt.output = false;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda2s = 1; opt.k = 30;
% 25 features
opt.lambda1s = 0.125; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 50 features
opt.lambda1s = 0.09; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 100 features
opt.lambda1s = 0.055; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));
% 200 features
opt.lambda1s = 0.03; smodel = DISMvis1(data, opt);
DISMauc(smodel{1}.idxSeleDISMtheta, data.gndFea, length(unique(data.testSet)));




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 8, Liver Metasis Dataset, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{4};

description = 'Exp-8-Liver-Metasis-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = [40,60,80,100];
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = logspace(0,-2,15);
opt.lambda2s =  [0.1, 0.5, 1,5,10, 20, 50, 100, 500]; 
% opt.lambda1s = 0.025:-0.002:0.005;
% opt.lambda2s = [0, 0.1, 0.5, 1:2:9, 20, 50, 100, 200, 500, 1000]; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);


[U, S, V]=SVD_CutOff(data.X,.95);
data.X = U'*data.X;
mSVM = accuracySVM(data); 
disp(sprintf('acc: %6.4f, std: %6.4f', mean(mSVM.accuCV), std(mSVM.accuCV)));



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 9, Breast Cancer Dataset, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{5};

description = 'Exp-9-Breast-Cancer-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = data.X(1:8141,:);
data.W = data.W(1:8141, 1:8141);
data.gndFea = data.gndFea(1:8141);
data.gName = data.gName(1:8141);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = 160;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = logspace(0,-2,10);
opt.lambda2s = 0;
% opt.lambda1s = 0.025:-0.002:0.005;
% opt.lambda2s = [0, 0.1, 0.5, 1:2:9, 20, 50, 100, 200, 500, 1000]; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 10, Maize Dataset, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load ppiData.mat
data = dataAll{8};

description = 'Exp-10-Maize-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);
% 
% data.X = data.X([1:19998,20000:end],:);
% data.W = data.W([1:19998,20000:end], [1:19998,20000:end]);
% data.gndFea = data.gndFea([1:19998,20000:end]);
% data.gName = data.gName([1:19998,20000:end]);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = [50,100];
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = logspace(0,-2,10);
% opt.lambda2s = 0;
% opt.lambda1s = 0.025:-0.002:0.005;
opt.lambda2s = [0, 0.1, 0.5, 1, 5, 10, 20, 50, 100, 200, 500]; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);
save([data.figName,'.mat'], 'smodel');




% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 11, Embryonic Dataset, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load ppiData.mat
data = dataAll{6};

description = 'Exp-11-Embry-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);
% 
% data.X = data.X([1:19998,20000:end],:);
% data.W = data.W([1:19998,20000:end], [1:19998,20000:end]);
% data.gndFea = data.gndFea([1:19998,20000:end]);
% data.gName = data.gName([1:19998,20000:end]);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = [20,30,50];
opt.nFeaUpd = .8; opt.verbose = 1;
% opt.lambda1s = logspace(0,-2,10);
opt.lambda1s = logspace(-3,-7,10);
opt.lambda2s = 0;
% opt.lambda1s = 0.025:-0.002:0.005;
% opt.lambda2s = [0, 0.5, 1, 5, 10,50]; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);
save([data.figName,'.mat'], 'smodel');












































