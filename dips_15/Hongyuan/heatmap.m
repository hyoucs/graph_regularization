clear all; close all; clc;
addpath('../');
load tempData.mat      % with lambda1 = 0.062
data = dataAll{3};

description = 'Exp-1-brain-DISM';
datasetName = strrep([description,'-',data.name,'-',datestr(now)], ' ', '-');
datasetName = strrep(datasetName, ':', '-');
datasetName = strrep(datasetName, '/', '');
% datasetName = strrep([description,'-',datestr(now)], ' ', '-');
mkdir([pwd, '/dataResult/', datasetName]);

data.X = Xnorm(data.X,1);
data.posetype = 'none';
data.figName = [pwd, '/dataResult/', datasetName];

opt.visfold = false; 
opt.visall = false;

opt.bLinear = 0; opt.alpha = .2; opt.k = 10;
opt.nFeaUpd = .8; opt.verbose = 1;

% opt.lambda1s = [0.0150:-0.005:0.005, 0.001];
opt.lambda1s = [0.5000, 0.3162, 0.3000, 0.2000, 0.1540, 0.1500, 0.1000, 0.0750,...
                 0.0500, 0.0365, 0.0178, 0.0100, 0.0087, 0.0042, 0.0021, 0.0010, 0];
opt.lambda2s = [0 0.1 0.5 1:2:9, 10:10:50, 70:20:150, 200:50:500, 100:500:2000, 3000];
% opt.lambda2s = logspace(log(0.0001)/log(10), log(0.1)/log(10), 10);

accMatrix = zeros(length(opt.lambda1s), length(opt.lambda2s));
smodel = DISMvis1(data, opt);
for i = 1:length(opt.lambda2s),
    for j = 1:length(opt.lambda1s),
        k = length(opt.lambda1s)*(i-1)+j;
        accMatrix(i,j) = mean(smodel{k}.accuCV);
        if accMatrix(i,j) < 0.2
            break;
        end
    end
end
save('heatmap-10.mat', 'accMatrix');


opt.lambda1s = [0.5000, 0.3162, 0.3000, 0.2000, 0.1540, 0.1500, 0.1000, 0.0750,...
                 0.0500, 0.0365, 0.0178, 0.0100, 0.0087, 0.0042, 0.0021, 0.0010, 0];
opt.lambda2s = [0 0.1 0.5 1:2:9, 10:10:50, 70:20:150, 200:50:500, 600:500:2000, 3000];
load heatmap-110.mat

accMatrix = accMatrix([1:25,27:30],:);
accMatrix = accMatrix(1:length(opt.lambda2s), 1:length(opt.lambda1s));

figure;
plot(opt.lambda1s, accMatrix(1,:), 'ro-'); hold on;
for i = 2:length(opt.lambda2s),
    plot(opt.lambda1s, accMatrix(i,:)); hold on;
end
title('knn=110')

[x, y] = meshgrid(opt.lambda1s, opt.lambda2s);
surf(x, y, accMatrix);
xlabel('lambda1'); ylabel('lambda2'); zlabel('accuracy');
title('accuracy surface vs parameter setting');
set(gca,'xscale','log','yscale','log');



figName = 'heatmap-knn-100';
saveas(gcf, [figName, '.pdf'], 'pdf');
saveas(gcf, [figName, '.fig'], 'fig');
saveas(gcf, [figName, '.png'], 'png');

