% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 5, Brain Dataset, curve and barplot, parameter tuning
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load tempData.mat      % with lambda1 = 0.062
data = dataAll{3};

description = 'Exp-5-brain-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'brain';
data.figName = [pwd, '/dataResult/', datasetName];

opt.visfold = false; opt.visall = false ;
opt.bLinear = 0; opt.alpha = .2; opt.k = 40;
opt.nFeaUpd = .8; opt.verbose = 1;
% opt.lambda1s = 0.029; opt.svm = false;
opt.lambda1s = [0.3:-0.01:0.22, 0.2:-0.005:0.01];
opt.lambda2s = 0.1;
opt.svm = true;
tStart = tic; m = DISMvis1(data, opt); toc(tStart);


dresAll = {};
res.acc = result(:,6)';         res.sd = result(:,7)';
res.size = result(:,4)';        res.culsize = result(:,5)';
res.lambda1 = result(:,1)';	 res.lambda2 = result(:,2)';    res.knn = result(:,3)';
res.legend = 'DISM w/o (k-40)';
dresAll{end+1} = res;

hf = figure;
han = {}; leg = {};
colorIndex = 1;
colorArray = [0 153 204; 0 0 255; 255 153 0; 255 0 0; ...
	153 0 153; 255 51 204; 153 51 0; 0 153 0; 0 0 0]./255;
lineArray = {'d-', '*-', '>-', '<-', 'd-', 's-'};
for i = 1:length(dresAll)
    testing_accuracy_size = dresAll{i}.size;
    testing_accuracy = dresAll{i}.acc;
	leg{end+1} = dresAll{i}.legend;
    han{end+1} = plot(ceil(testing_accuracy_size), testing_accuracy,...
        char(lineArray(mod(colorIndex-1,length(lineArray))+1)),...
        'color', colorArray(colorIndex,:), 'LineWidth', 1.5);
    colorIndex = colorIndex + 1;
    hold on;
end
han{end+1} = plot(1:300, ones(1,300).*0.5530, '-',...
    'color', colorArray(colorIndex,:), 'LineWidth', 1.5);
axis([0, 120, 0.5, 1]);
h = []; for i = 1:length(han), h = [h, han{i}]; end
hl = legend(h, leg{1},leg{2},'SVM', 2);
xlabel('number of features'); ylabel('classification accuracy');
set(hf, 'Position', [100 100 600 500]);



% PCA-SVM
[U, S, V]=SVD_CutOff(data.X,.95);
data.X = U'*data.X;
mSVM = accuracySVM(data);
disp(sprintf('acc: %6.4f, std: %6.4f', mean(mSVM.accuCV), std(mSVM.accuCV)));

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Step 6, Brain Dataset, common selected features
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
addpath('../');
load tempData.mat      % with lambda1 = 0.062
data = dataAll{3};

description = 'Exp-6-brain-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'brain';
data.figName = [pwd, '/dataResult/', datasetName];

opt.visfold = false; opt.visall = false;
opt.bLinear = 0; opt.alpha = .2; 
opt.nFeaUpd = .8; opt.verbose = 1;

% DISM
opt.lambda1s = 0.0350; opt.lambda2s = 0.1000; opt.k = 40;
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);

% DISM w/o
opt.lambda1s = 0.0300; opt.lambda2s = 0; opt.k = 40;
tStart = tic; smodel = DISMvis1(data, opt); toc(tStart);

% compute common features
idxSele = [];
for i = 1:length(unique(data.testSet))
	if isempty(idxSele)
		idxSele = smodel{1}.idxSeleDISMfold{i};
	end
	idxSele = intersect(idxSele, smodel{1}.idxSeleDISMfold{i});
	disp(['Fold #',num2str(i),': ',num2str(length(smodel{1}.idxSeleDISMfold{i})),' features']);
end
disp(['Common Features: ',num2str(length(idxSele))]);

% transfer selected nodes to selected edges
name = data.cellName;
edgeMat = triu(gen_map(length(name)));
for i = 1:length(idxSele)
    [r, c] = find(edgeMat == idxSele(i));
    fprintf('%4i: %s -- %s \n',idxSele(i), name{r}, name{c});
    for j = 1:length(unique(data.testSet))
        [~, I] = sort(abs(smodel{1}.idxSeleDISMtheta{j}), 'descend');
        rank = find(I == idxSele(i));
        fprintf('%d, ', rank);
    end
    fprintf('\n');
end




