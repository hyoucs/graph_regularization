% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment 1: fix lambda1 = 0.07 and k = 50, change lambda2, increase the
% number of selected features and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
% load lambda2effect.mat      % with lambda1 = 0.07
load noised-straight-pose.mat      % with lambda1 = 0.062

% create a folder for a dataset
description = 'Experiment 1';
% datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
% datasetName = strrep([datasetName], ':', '-');
% datasetName = strrep([datasetName], '/', '');
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.lambda1s = 0.062;
opt.lambda2s = [0 0.1 0.5, 1:10, 20:20:100, 200:200:1000];
opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
% supLE = supLE(data,opt);
tStart = tic; m = DISMvisual(data, opt); toc(tStart) % m.accuCV

mSVM = accuracySVM(data);
disp(strcat('avg svm accuracy:',num2str(mean(mSVM.accuCV)), '(' ,num2str(std(mSVM.accuCV)), ')'));



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment 2: fix k = 50, change lambda1 and lambda2 at the same time,
% fix the number of selected features around 20 and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
% load lambda2effect.mat      % with lambda1 = 0.07
load noised-straight-pose.mat% with lambda1 = 0.062

% create folder for a dataset
description = 'Experiment 2';
% datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
% datasetName = strrep([datasetName], ':', '-');
% datasetName = strrep([datasetName], '/', '');
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.figName = [pwd, '/dataResult/', datasetName];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;

opt.lambda1s = 0.062; opt.lambda2s = 0;
tStart = tic; data.posetype = 'straight';
m = DISMvisual(data, opt); data.posetype = 'none'; toc(tStart)
feaNum = m.feaNum;

lambda2Array = [0.1 0.5, 1:10, 20:20:100, 200:200:1000];
for i = 1:length(lambda2Array)
    tuning = true;
    data.posetype = 'none';
    opt.lambda1s = 0.062;
    while tuning
        opt.lambda2s = lambda2Array(i);
        tStart = tic; m = DISMvisual(data, opt); toc(tStart)
        feaNumCurr = m.feaNum;
        if (feaNumCurr - feaNum > 1)
            opt.lambda1s = opt.lambda1s + 0.005;
        else
            tuning = false;
            data.posetype = 'straight';
            m = DISMvisual(data, opt);
            disp('visualization figure saved'); disp(' ');
        end
    end
end

mSVM = accuracySVM(data);
disp(strcat('avg svm accuracy:',num2str(mean(mSVM.accuCV)), '(' ,num2str(std(mSVM.accuCV)), ')'));


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment 3: fix k = 50, fix lambda2 = 0, and change lambda1,
% and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
load morenoise-19-9.mat

description = 'Experiment 3';
% datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
% datasetName = strrep([datasetName], ':', '-');
% datasetName = strrep([datasetName], '/', '');
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = [0.3:-0.02:0.07, 0.06:-0.01:0.01];
opt.lambda2s = [0 0.5];
tStart = tic; m = DISMvisual(data, opt); toc(tStart);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment 4: fix lambda1 = 0.018, fix lambda2 = 0, and change k,
% and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
load lambda2effect.mat

description = 'Experiment 4';
% datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
% datasetName = strrep([datasetName], ':', '-');
% datasetName = strrep([datasetName], '/', '');
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.bLinear = 0; opt.alpha = .2;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = 0.018;
opt.lambda2s = 0;

kArray = [30:5:50, 60:10:150, 170:30:300];
for i = 1:length(kArray)
    opt.k = kArray(i);
    tStart = tic; m = DISMvisual(data, opt); toc(tStart);
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment 5: fix lambda1 = 0.018, fix k = 50, and change lambda2,
% and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
load lambda2effect.mat

description = 'Experiment 5';
% datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
% datasetName = strrep([datasetName], ':', '-');
% datasetName = strrep([datasetName], '/', '');
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = 0.018;
opt.lambda2s = [0:0.01:1, 2:10, 20:20:100, 200:200:1000];
tStart = tic; m = DISMvisual(data, opt); toc(tStart);




% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment step 1: fix k = 50, change lambda1 and lambda2 at the same time,
% fix the number of selected features around 20 and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
load tempData.mat      % with lambda1 = 0.062
data = dataAll{1};

% create folder for a dataset
description = 'step-1-DISM';
% datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
% datasetName = strrep([datasetName], ':', '-');
% datasetName = strrep([datasetName], '/', '');
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.figName = [pwd, '/dataResult/', datasetName];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;

opt.lambda1s = 0.02; opt.lambda2s = 0.1;
tStart = tic; data.posetype = 'straight';
m = DISMvisual(data, opt); data.posetype = 'none'; toc(tStart)


feaNum = 100;
lambda2Array = [0.1 0.5, 1:10, 20:20:100, 200:200:1000];
for i = 1:length(lambda2Array)
    tuning = true;
    data.posetype = 'none';
    opt.lambda1s = 0.02;
    while tuning
        opt.lambda2s = lambda2Array(i);
        tStart = tic; m = DISMvisual(data, opt); toc(tStart)
        feaNumCurr = m.feaNum;
        if (feaNumCurr - feaNum > 1)
            opt.lambda1s = opt.lambda1s + 0.002;
        else
            tuning = false;
            data.posetype = 'straight';
            m = DISMvisual(data, opt);
            disp('visualization figure saved'); disp(' ');
        end
    end
end

mSVM = accuracySVM(data);
disp(strcat('avg svm accuracy:',num2str(mean(mSVM.accuCV)), '(' ,num2str(std(mSVM.accuCV)), ')'));



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
% Experiment 3 step 1: fix k = 50, fix lambda2 = 0, and change lambda1,
% and visualize selected subnetworks
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

clear all; close all; clc;
addpath('../');
% rmpath(genpath('./libsvm'));
load tempData.mat      % with lambda1 = 0.062
data = dataAll{1};

description = 'Exp-1-straight-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s = [0.175, 0.150, 0.105];
% opt.lambda1s = [0.3:-0.01:0.22, 0.2:-0.005:0.01];
opt.lambda2s = 350;
tStart = tic; DISMvis1(data, opt); toc(tStart);


% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
res.acc = result(:,6)';         res.sd = result(:,7)';
res.size = result(:,4)';        res.culsize = result(:,5)';
res.lambda1 = result(:,1)';	 res.lambda2 = result(:,2)';    res.knn = result(:,3)';
res.legend = 'DISM w/o (k-80)';    
res.name = 'DISM w/o';
dresAll = {};

if exist('all-poses-noise-dism.mat','file') == 0
	dresAll = {};
	save('all-poses-noise-dism.mat', 'dresAll');
end
load all-poses-noise-dism.mat
dresAll{end+1} = res;
save('all-poses-noise-dism.mat', 'dresAll');

if exist('brain-dism.mat','file') == 0
	dresAll = {};
	save('brain-dism.mat', 'dresAll');
end
load brain-dism.mat
dresAll{end+1} = res;
save('brain-dism.mat', 'dresAll');















clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'Exp-2-CMU-all-DISM';
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

feaNum = 75;
lambda2Array = [0.1 0.5, 1:10, 20:20:100, 200:200:1000];
for i = 1:length(lambda2Array)
    tuning = true;
    data.posetype = 'none';
    opt.lambda1s = 0.025;
    while tuning
        opt.lambda2s = lambda2Array(i);
        tStart = tic; m = DISMvis1(data, opt); toc(tStart)
        feaNumCurr = m{1}.feaNum;
        if (feaNumCurr - feaNum > 1)
            opt.lambda1s = opt.lambda1s + 0.005;
        else
            tuning = false;
            data.posetype = 'straight';
            m = DISMvis1(data, opt);
            disp('visualization figure saved'); disp(' ');
        end
    end
end


