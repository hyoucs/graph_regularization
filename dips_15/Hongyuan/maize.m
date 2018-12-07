function maize

addpath('../');
load ppiData.mat
data = dataAll{8};

description = 'DBL-Exp-10-Maize-DISM';
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
opt.lambda2s = 0;
% opt.lambda1s = 0.025:-0.002:0.005;
% opt.lambda2s = [0, 0.1, 0.5, 1, 5, 10, 20, 50, 100, 200, 500]; 
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);
save([data.figName,'.mat'], 'smodel');


end