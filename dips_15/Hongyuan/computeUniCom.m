% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 2, Brain dataset
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% DISMw/o
% acc    0.6470    0.6930    0.7460    0.7690    0.7570
% std    0.0310    0.0350    0.0410    0.0440    0.0580
% #nodes    9.2000   31.2000   48.2000   69.0000   87.0000

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{3};

description = 'Exp-2-Brain-DISMwo';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
% mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = 40;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s =  [0.1100,0.0650,0.0460,0.0300,0.0200];
opt.lambda2s = 0;
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);

for i = 1:length(smodel)
	uni = length(unique(smodel{i}.idxSeleDISM));
	com = length(unique(smodel{i}.idxSeleDISMcom));
    fprintf('%5.4f,%5.4f    ', com, uni);
end
fprintf('\n');

load 'DISMwo.mat'
smodel.model = smodel;
smodel.dataset = datasetName;
modelAll{end+1} = smodel;
save('DISMwo.mat', 'modelAll');



% DISM
% acc    0.6590    0.7290    0.7460    0.7810    0.7980
% std    0.0230    0.0730    0.0370    0.0510    0.0340
% #nodes    9.6000   30.6000   51.8000   70.8000   92.4000

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{3};

description = 'Exp-2-Brain-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
% mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = 40;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s =  [0.1150,0.0750,0.0550,0.0450,0.0350];
opt.lambda2s = 0.1;
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);

for i = 1:length(smodel)
	uni = length(unique(smodel{i}.idxSeleDISM));
	com = length(unique(smodel{i}.idxSeleDISMcom));
    fprintf('%5.4f,%5.4f    ', com, uni);
end
fprintf('\n');

load 'DISM.mat'
smodel.model = smodel;
smodel.dataset = datasetName;
modelAll{end+1} = smodel;
save('DISM.mat', 'modelAll');