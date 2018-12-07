% Compute precision and recall

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 1, CMU straight-poses, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% DISMw/o
% acc    0.7820    0.7820    0.7110    0.6930    0.6670
% std   0.0750    0.0590    0.0400    0.0780    0.0930
% #nodes    9.8000   28.8000   50.2000   73.2000   93.0000

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'Exp-1-CMU-straight-DISMwo';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
% mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s =  [0.1150, 0.0600, 0.0400, 0.0250, 0.0150];
opt.lambda2s = 0;
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);

prec = []; rec = [];
for i = 1:length(smodel)
	prec = [prec mean(smodel{i}.idxSeleDISMprec)];
	rec = [rec mean(smodel{i}.idxSeleDISMrec)];
    fprintf('%5.4f,%5.4f    ', prec(i), rec(i));
end
fprintf('\n');

load 'DISMwo.mat'
smodel.model = smodel;
smodel.dataset = datasetName;
modelAll{end+1} = smodel;
save('DISMwo.mat', 'modelAll');



% DISM
% acc    0.7250    0.8070    0.8140    0.7820    0.7880
% std    0.0340    0.0770    0.0740    0.0710    0.1140
% #nodes   10.8000   29.0000   48.8000   68.2000   84.6000

clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'Exp-1-CMU-straight-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
% mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'straight';
data.figName = [pwd, '/dataResult/', datasetName];

opt.svm = true; opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = .8; opt.verbose = 1;
opt.lambda1s =  [0.2000, 0.1400, 0.1050, 0.0900, 0.0800];
opt.lambda2s = 350;
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);

prec = []; rec = [];
for i = 1:length(smodel)
	prec = [prec mean(smodel{i}.idxSeleDISMprec)];
	rec = [rec mean(smodel{i}.idxSeleDISMrec)];
    fprintf('%5.4f,%5.4f    ', prec(i), rec(i));
end
fprintf('\n');

load 'DISM.mat'
smodel.model = smodel;
smodel.dataset = datasetName;
modelAll{end+1} = smodel;
save('DISM.mat', 'modelAll');


% visualize with gndFea
clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'DISM-CMU-prec-rec-check';
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'gndsep';
data.figName = [pwd, '/dataResult/', datasetName];

opt.visfold = true;
opt.visall = true;

load DISM.mat
idxModelAll = 1;
for idx_model = 1:length(modelAll{idxModelAll}.model)
	m = modelAll{idxModelAll}.model{idx_model};
    for idx_fold = 1:length(unique(data.testSet))
        if isfield(opt, 'visfold') && opt.visfold
            figName = sprintf('%s/fold_%d_%6.3f_%6.3f_%d_fea_%6.3f_%6.3f', ...
                data.figName, idx_fold, m.lambda1, m.lambda2, m.k, ...
                length(m.idxSeleDISMfold{idx_fold}), m.acc(idx_fold));
            visual(data, m.idxSeleDISMfold{idx_fold}, figName);
        end
    end
    if isfield(opt, 'visall') && opt.visall
        figName = sprintf('%s/%6.3f_%6.3f_%d_avgFea_%6.3f_culFea_%d_%6.3f', ...
                data.figName, m.lambda1, m.lambda2, m.k, ...
                mean(m.idxSeleDISMsize), length(m.idxSeleDISM), m.accuCV);
        visual(data, m.idxSeleDISM, figName);
    end
end




% Compute precision and recall with union features

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 1, CMU straight-poses, curve, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


% DISM
clear all; close all; clc;
addpath('../');
load tempData.mat
data = dataAll{1};

description = 'DISM-CMU-prec-rec-check';
datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'gnd';
data.figName = [pwd, '/dataResult/', datasetName];

opt.visfold = true;
opt.visall = true;

load DISM.mat
idxModelAll = 1;
for idx_model = 1:length(modelAll{idxModelAll}.model)
	m = modelAll{idxModelAll}.model{idx_model};
    for idx_fold = 1:length(unique(data.testSet))
        if isfield(opt, 'visfold') && opt.visfold
            figName = sprintf('%s/fold_%d_%6.3f_%6.3f_%d_fea_%6.3f_%6.3f', ...
                data.figName, idx_fold, m.lambda1, m.lambda2, m.k, ...
                length(m.idxSeleDISMfold{idx_fold}), m.acc(idx_fold));
            visual(data, m.idxSeleDISMfold{idx_fold}, figName);
        end
    end
    if isfield(opt, 'visall') && opt.visall
        figName = sprintf('%s/%6.3f_%6.3f_%d_avgFea_%6.3f_culFea_%d_%6.3f', ...
                data.figName, m.lambda1, m.lambda2, m.k, ...
                mean(m.idxSeleDISMsize), length(m.idxSeleDISM), m.accuCV);
        visual(data, m.idxSeleDISM, figName);
    end
end

