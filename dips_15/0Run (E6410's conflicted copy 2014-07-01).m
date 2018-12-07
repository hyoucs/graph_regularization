

% -- feaSubnet with brain data
%%
clear all;close all; clc;
% %%
% load 'data/3Gauss150Dim'; whos
% load 'data/2Gauss150Dim'; whos
% load 'data/ppiData'; whos

load data/brain; whos
nClass = length(unique(data.gnd));
for iClass = 3:nClass
    load data/brain; whos
%     % remove each class from data
%     clsIdx = find(data.gnd==iClass);
%     data.X(:,clsIdx) =[];
%     data.gnd(:,clsIdx) =[];
%     data.testSet(:,clsIdx) =[];
    
    % merging classes
    cls1 = (data.gnd == iClass);
    cls2 = ~cls1;
    data.gnd(cls1) = 1;
    data.gnd(cls2) = 2;
    

    data.X = Xnorm(data.X,1);
%     opt.lambda1s = [logspace(2, 0, 7) 0]/1000;
    opt.lambda1s = [0.01];
    opt.lambda2s = [.5];

    opt.bLinear = 0; opt.alpha = .1; opt.k = 50; 
    opt.nFeaUpd = .5; opt.verbose = 1;
    % supLE = supLE(data,opt);
    tStart = tic;m = feaSubnet(data,opt); toc(tStart) % m.accuCV
    mSVM = accuracySVM(data); 
    disp(strcat('avg svm accuracy:',num2str(mean(mSVM.accuCV)), '(' ,num2str(std(mSVM.accuCV)), ')'));
    
end


option.topFea=find(m.theta~=0);
subGraphDisp(data.W,option);

%-- disp edge name here



% --- LDA and supLE
%%
clear all;close all; clc;
% %%
% load 'data/3Gauss150Dim'; whos
% load 'data/2Gauss150Dim'; whos
load data/brain; whos
data.X = Xnorm(data.X,1);
% m = accuracySVM(data);  %accuCV: [0.5902 0.5161 0.5738 0.6066 0.4754]

tic
[Y,model] = lda(data.X',data.gnd,length(unique(data.gnd))-1);
toc

tic
opt.bLinear = 1; opt.alpha = .1; opt.k = 50; opt.lowerBound = 80;
opt.nFeaUpd = .5; opt.verbose = 1;
supLE = supLE(data,opt);
toc



%% ------- svm test - done
clear all;close all; clc;
load 'data/3Gauss150Dim'; whos
m = accuracySVM(data)
m.accuCV


%% ----------  LDA

clear all;close all; clc;
load 'data/prostateData'; whos
% load 'data/3Gauss150Dim'; whos
% load 'data/2Gauss150Dim'; whos
% load data/brain; whos
tic
[Y,model] = lda(data.X',data.gnd,length(unique(data.gnd))-1);
toc




%% test eVecs returned by LDA
% Concl: 2 top coeff in 2 eVecs are most important
% yet, accuracy is high only when at least top 6000 are kept!!!


clear all;close all; clc;
load data/brain;
load data/brainRes01; whos
m = LDAmodel;
nTopFea = 6000;

for i = 1:size(m.eVecs,2)
    [dmp idx] = sort(abs(m.eVecs(:,i)),'descend');
    m.eVecs(idx(nTopFea+1:end),i)=0; % remove less important coeff
    %     m.eVecs(idx(1:nTopFea),i)=0; % remove most important coeff
end

figure;
subplot(1,2,1); stem(abs(m.eVecs(:,1)));
subplot(1,2,2); stem(abs(m.eVecs(:,2)));

mappedX = m.eVecs'*data.X;

figure;
subplot(1,2,1); PlotX(data.X,data.gnd,'','Original space',''); grid on;
subplot(1,2,2); PlotX(mappedX,data.gnd,'','Projected space','');

dataProj.X = mappedX;
dataProj.gnd = data.gnd;
dataProj.testSet = data.testSet;

m = accuracySVM(dataProj);
m.accuCV

%-- end test with LDA's eigenVectors


%%


load data/brain; whos
[nFea nSmp] = size(data.X);
data.X = Xnorm(data.X);

% get y from supLE
opt.bLinear=0; opt.alpha=.1;opt.k=50;
m1 = supLE(data,opt);


figure;
subplot(1,2,1);PlotX(m1.Y,data.gnd,'','','');
subplot(1,2,2);stem(abs(m1.eVecs(:,1)));

%%
data.y = Xnorm(m1.Y(1,:));
% data.y=data.gnd;

wLS = data.X'\m1.Y(1,:)'; denom = sum(abs(wLS'));
lambda1s = [logspace(3, 0, 50) 0];
lambda1s = lambda1s/100;

w=zeros(nFea, length(lambda1s));
tic;
options.lambda2=10.9;
for i=1:length(lambda1s)
    options.lambda1=lambda1s(i);options.verbose=0;
    m = regCoordAscent(data, options); w(:,i) = m.a;  %-- my own
    %     w(:,i) = eNet(data,options); %-- elasticNet with network
end
fprintf('Overall Elapsed time: %8.2f \n\n',toc);

w = w';
s = sum(abs(w),2)/denom;
figure;plot(s, w, '-o');
title(strcat('lambda2=',num2str(options.lambda2)));
% legend(data.fName, 'location', 'northwest');
% figure;imshow(w'*1000,'InitialMagnification',2000);
% w(15:end,:)'
%%





%% DONE with sequential approach: (1) (nonlinear)supLE; (2) (linear)regCoordAscent
clear all;close all; clc;
% %%
load 'data/prostateData'; whos
% load 'data/2Gauss150Dim'; whos
% load 'data/toyExamp';

data.X = Xnorm(data.X);
[nFea nSmp] = size(data.X);

% get y from supLE
opt.bLinear=0; opt.alpha=.1;opt.k=50;
m1 = supLE(data,opt);
% figure;
% subplot(1,2,1);PlotX(m1.Y,data.gnd,'','','');
% subplot(1,2,2);stem(abs(m1.eVecs(:,1)));

data.y = Xnorm(m1.Y(1,:));
% data.y=data.gnd;

wLS = data.X'\m1.Y(1,:)'; denom = sum(abs(wLS'));
lambda1s = [logspace(3, 0, 50) 0];
lambda1s = lambda1s/100;

w=zeros(nFea, length(lambda1s));
tic;
options.lambda2=0.5;
for i=1:length(lambda1s)
    options.lambda1=lambda1s(i);options.verbose=0;
    m = regCoordAscent(data, options); w(:,i) = m.a;  %-- my own
    %     w(:,i) = eNet(data,options); %-- elasticNet with network
end
fprintf('Overall Elapsed time: %8.2f \n\n',toc);

w = w';
s = sum(abs(w),2)/denom;
figure;plot(s, w, '-o');
title(strcat('lambda2=',num2str(options.lambda2)));
% legend(data.fName, 'location', 'northwest');
figure;imshow(w'*1000,'InitialMagnification',2000);
% w(15:end,:)'
%%


%%
%-------------------------------------------------------------------------%
%-----------------  coordinateAscentLDA or (pLDA)  -----------------------%
%-- max v'Bv is meaningfull given its definition
%-------------------------------------------------------------------------%
%-- 2Gauss150Dim: pLDA (LDA ->coorodinateAscent) with lambda1s/1000 and check w(:,24:27)


clear all;close all; clc;
load 'data/prostateData'; whos
% load 'data/2G20Dim3GT'; whos
% load 'data/2Gauss150Dim'; whos
% load 'data/3Gauss150Dim'; whos
%%
% load 'data/toyExamp';
% data.X(6,:) = data.X(1,:);
% data.X(7,:) = data.X(1,:);


data.X = Xnorm(data.X);
[nFea nSmp] = size(data.X);

% %-- LDA for B,W, initial theta
% [Y,m0] = lda(data.X',data.gnd);
% figure;
% subplot(1,2,1);PlotX(Y',data.gnd,'','','');
% subplot(1,2,2);stem(abs(m0.eVecs(:,1)));
% data.B = m0.B; data.W = m0.W;


%-- get y from supLE
opt.bLinear=0; opt.alpha=.1;opt.k=50;
m1 = supLE(data,opt);

% figure;
% subplot(1,2,1);PlotX(m1.Y,data.gnd,'','','');
% subplot(1,2,2);stem(abs(m1.eVecs(:,1)));

% data.y = m1.Y(1,:);
data.y = Xnorm(m1.Y(1,:));
% data.y=data.gnd;

wLS = data.X'\m1.Y(1,:)'; denom = sum(abs(wLS'));
lambda1s = [logspace(3, 0, 50) 0];
lambda1s = lambda1s/100;

w=zeros(nFea, length(lambda1s));

options.lambda2=0.5;
% for i=20:20%length(lambda1s)
for i=15:length(lambda1s)
    options.lambda1=lambda1s(i);
    m = regCoordAscent(data, options); w(:,i) = m.a;  %-- my own
    %     w(:,i) = eNet(data,options); %-- elasticNet with network
end

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;plot(s1, w1, '-o');
title(strcat('lambda2=',num2str(options.lambda2)));
legend(data.fName, 'location', 'northwest');

w(:,15:end)


%%

%-- linear regression
wLS = data.X'\data.gnd'; denom = sum(abs(wLS'));

lambda2 = 1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda1s = lambda1s/100;
lambda2s = lambda2*ones(1,length(lambda1s));

% lambda1s = [22.29 22];


w=zeros(nFea, length(lambda1s));
options.theta = m1.eVecs(:,1);
options.theta = data.X'\m1.eVecs(:,1);
%-- coordinateAscent
for i=1:length(lambda1s)
    options.lambda1=lambda1s(i);
    m = coordinateAscent(data,options);
    w(:,i) = m.v;
end

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;
subplot(1,2,1);plot(s1, w1, '-o');
subplot(1,2,2);imshow(w,'InitialMagnification',2000);





%%

%-------------------------------------------------------------------------%
%--------------------------RegSubnetMajorization--------------------------%
%-------------------------------------------------------------------------%



clear all;close all; clc;
load 'data/prostateData'; whos

% load 'data/2Gauss150Dim'; whos



[nFea nSmp] = size(data.X);
data.W = rand(nFea,nFea);

wLS = data.X'\data.gnd'; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda2s = lambda2*ones(1,length(lambda1s));


options.alpha=0;
options.k=15;
options.dDim=1;

w=zeros(nFea, length(lambda1s));
for i=1:length(lambda1s)
    options.beta=lambda1s(i);
    m = RegSubnetMajorization(data, options, lambda1s(i)); %-- my adapted shooting algo
    w(:,i) = m.u;
end

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;plot(s1, w1, '-o');
title('Elastic net on prostate data');
legend(data.feaName, 'location', 'northwest');
% set(gca,'ylim',[-.5 1])
% set(gca,'ylim',[-2 3]);set(gca,'xlim',[0 1]);
%%


%%
%-------------------------------------------------------------------------%
%-----------------  penalizedLDA & SparseLDA  ----------------------------%
%-------------------------------------------------------------------------%
%-- 2Gauss150Dim: penalizedLDA runs perferctly with lambda1 >= 0.02
%-- sparseLDA: need to look back!

clear all;close all; clc;
% load 'data/prostateData'; whos
% load 'data/2Gauss150Dim'; whos
load 'data/3Gauss150Dim'; whos
% load 'data/toyExamp'; whos


[nFea nSmp] = size(data.X);
data.W = rand(nFea,nFea);

wLS = data.X'\data.gnd'; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda1s = lambda1s/1000;
lambda2s = lambda2*ones(1,length(lambda1s));


%-- Witten penalizedLDA
for i=1:length(lambda1s)
    options.lambda1=lambda1s(i);
    m = penalizedLDA(data, options, lambda1s(i));
    w(:,i) = m.v;
end
%-- End Witten penalizedLDA

%
% %-- Clemmensen SparseLDA
% uniClass = unique(data.gnd);
% nClass = length(uniClass);
% Y = zeros(nSmp,nClass);
% for idxCls=1:nClass
%     idxSmp = find(data.gnd==uniClass(idxCls));
%     Y(idxSmp,idxCls) = 1;
% end
% X=data.X';
% for i=1:length(lambda1s)
%   [beta os] = slda(X,Y, lambda1s(i), 0, 1, 1000);
%   w(:,i) = beta;
% end
% %-- End Clemmensen SparseLDA

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;subplot(1,2,1);plot(s1, w1, '-o');
subplot(1,2,2);imshow(w,'InitialMagnification',2000);
% set(gca,'ylim',[-.5 1])
% set(gca,'ylim',[-2 3]);set(gca,'xlim',[0 1]);







%%
%-------------------------------------------------------------------------%
%---------------------------- supLE DONE! ---------------------------%
%-------------------------------------------------------------------------%
%-- works well on syn with both linear/non-linear cases

% 2Gauss150Dim: sensitive to k setting, see separation by k=5 and k=20
% 3Gauss150Dim: likewise, k should be at least 45 (larger, better)
%               for non-linear method, even much better than linear one with k>50

clear all;close all; clc;
load 'data/2Gauss150Dim';
% load 'data/3Gauss150Dim';

options.alpha=0; options.k=50; options.bLinear=1;
m = supLE(data,options);
Y = m.Y; [data.gnd' Y']

figure;
subplot(1,2,1); PlotX(data.X,data.gnd,'','',''); grid on;
subplot(1,2,2); PlotX(Y,data.gnd,'','',''); grid on;
%%
%-- prostateData: supLE works well
% 2 steps: 1. learn y as a non-linearly classifying boundary
%           (or learn y as linear proj of X on eVec from LSDA)
%          2. approx y by eNet (shooting algo)
% Results:
%   + both linear and nonlinear cases resemble L1Shooting since
%       y resembles the class labels
%   + for NONlinear case: need to set k large (say, 30)!

% clear all;close all; clc;
load 'data/prostateData'; whos
% data.gnd=sign(data.gnd);
% figure;PlotX(data.X,data.gnd,'','',''); grid on;
options.alpha=0; options.k=55; options.bLinear=0;
model = supLE(data,options);


Y=Xnorm(model.Y(1,:),1);
[nFea nSmp] = size(data.X);
data.W = rand(nFea,nFea);
wLS = data.X'\Y'; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda2s = lambda2*ones(1,length(lambda1s));
w=zeros(nFea, length(lambda1s));
for i=1:length(lambda1s)
    options.beta=lambda1s(i);
    w(:,i) = L1Shooting(data.X', Y', lambda1s(i)); %-- ori.Shooting
end
w1 = w';
s1 = sum(abs(w1),2)/denom;


figure;
subplot(1,2,1);plot(Y(data.gnd<0),'ro'); hold on; plot(Y(data.gnd>0),'b*');
subplot(1,2,2); plot(s1, w1, '-o');
title('Elastic net on prostate data');
legend(data.fName(1:8), 'location', 'northwest');
set(gca,'xlim',[0 1]);% set(gca,'ylim',[-.5 1])

%%
%-------------------------------------------------------------------------%
%-- compared to this L1shooting algorithm
% clear all;close all; clc;
load 'data/prostateData'; whos

[nFea nSmp] = size(data.X);
data.W = rand(nFea,nFea);

wLS = data.X'\data.gnd'; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda2s = lambda2*ones(1,length(lambda1s));


options.alpha=0;
options.k=10;
options.dDim=1;

w=zeros(nFea, length(lambda1s));
for i=1:length(lambda1s)
    options.beta=lambda1s(i);
    w(:,i) = L1Shooting(data.X', data.gnd', lambda1s(i)); %-- my adapted shooting algo
end

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;plot(s1, w1, '-o');
title('Elastic net on prostate data');
legend(data.fName(1:8), 'location', 'northwest');
set(gca,'ylim',[-2 3]);set(gca,'xlim',[0 1]);

%%





%-------------------------------------------------------------------------%







%%
%-------------------------------------------------------------------------%
%---------------------------- elasticNetDemo -----------------------------%
%-------------------------------------------------------------------------%
% elasticNetDemo;

clear all;close all; clc;
load 'data/prostate'
y = (y>mean(y)); %-- convert to classification problem
X = Xnorm(X')';
X = mkUnitNorm(X);
y = Xnorm(y')';
[n p] = size(X);
wLS = X\y; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda2s = lambda2*ones(1,length(lambda1s));

xbar = mean(X);
XtrainC = X - repmat(xbar,size(X,1),1);
ybar = mean(y);
ytrainC = y-ybar;

for i=1:length(lambda1s)
    model = elasticNet(XtrainC, ytrainC, lambda1s(i), lambda2s(i));
    w0 = ybar - xbar*model;
    w(:,i) = [w0; model];
end


w1 = w(2:end,:)'; % skip offset , w1(iter, var)
s1 = sum(abs(w1),2)/denom;
figure;clf
plot(s1, w1, '-o');
%plot(log(lambda1s), w1, '-o');
title('Elastic net on prostate data')
legend(names(1:8), 'location', 'northwest')
set(gca,'ylim',[-2 3])
xlabel(sprintf('shrinkage factor s(%s)', '\lambda_1'))





%%
%-------------------------------------------------------------------------%
%---------------------------- L1Shooting -----------------------------%
%-------------------------------------------------------------------------%



clear all;close all; clc;
load 'data/prostate'
Y = (Y>mean(Y)); %-- convert to classification problem
X = Xnorm(X')';
X = mkUnitNorm(X);
Y = Xnorm(Y')';
[n p] = size(X);
wLS = X\Y; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda2s = lambda2*ones(1,length(lambda1s));

xbar = mean(X);
XtrainC = X - repmat(xbar,size(X,1),1);
ybar = mean(Y);
ytrainC = Y-ybar;

for i=1:length(lambda1s)
    model = L1Shooting(XtrainC, ytrainC, lambda1s(i));
    w(:,i) = model;
end

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;plot(s1, w1, '-o');
title('Elastic net on prostate data');
legend(names(1:8), 'location', 'northwest');
set(gca,'xlim',[0 1.1]);set(gca,'ylim',[-2 3]);


%%
%-- test on penalizedLDA by Witten and regSubnet
%-- 2Gauss150Dim: penalizedLDA runs perfectly with lambda1 >= 0.02


clear all;close all; clc;
load 'data/2Gauss150Dim'; whos

[nFea nSmp] = size(data.X);
data.W = rand(nFea,nFea);

wLS = data.X'\data.gnd'; denom = sum(abs(wLS'));
lambda2=1000;
lambda1s = [logspace(3, 0, 50) 0];
lambda1s = lambda1s/100;
lambda2s = lambda2*ones(1,length(lambda1s));


%-- Witten penalizedLDA
for i=1:length(lambda1s)
    options.lambda1=lambda1s(i); options.k=20;
    m = penalizedLDA(data, options, lambda1s(i));
    %     m = regSubnet(data, options, lambda1s(i));
    w(:,i) = m.v;
end

w1 = w';
s1 = sum(abs(w1),2)/denom;
figure;plot(s1, w1, '-o');
title('Elastic net on prostate data');
legend(data.fName(1:8), 'location', 'northwest');

w(1:30,24:28)


%%
% %-------------------------------------------------------------------------%
% %---------------------------- toy example -----------------------------%
% %-------------------------------------------------------------------------%
% %-- scenario setup:
% %-- nodes 1,2,3 are relevants, except 1 val in node 1 is not consistent
% %-- They thus should all selected to form subnetwork!
% %-- nodes 6, 7 are completely relevant yet no edge connecting, so should

% not be selected
data.X = [0 0 1 1 0; 1 0 1 1 0; 0 1 0 0 1; 1 1 0 1 0; 0 0 0 1 1; 1 0 1 1 0; 1 0 1 1 0];
data.gnd = [1 0 1 1 0];
data.W = [0 1 1 0 0 0 0; 1 0 0 1 0 0 0; 1 0 0 1 0 0 0; 0 1 1 0 1 0 0; 0 0 0 1 0 1 1; 0 0 0 0 1 0 0; 0 0 0 0 1 0 0];
data.feaName = {'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7'};
data.gndFea = [1 1 1 0 0 0 0];
save('data/toyExamp.mat','data');

%-- plot network
W=triu(data.W);
idxFea=[1:size(W,1)];
% feaName= strread(num2str(idxFea),'%s');
feaName= data.fName;


% load data/toyExamp
% W = data.W
gObj = biograph(W,feaName,'ShowArrows','off') % get(bg2.nodes,'ID')
gObj = view(gObj);

%-- color GTnodes
gtNodeIdx= find(data.gndFea==1); %-- find gtNodes in top 100 nodes
set(gObj.nodes(gtNodeIdx),'Color',[1 0 0],'size',[40 30]);
dolayout(gObj);


%% confusion mat
close all; clear all; clc;
load fisheriris
numObs = length(species);
p = randperm(numObs);
% swap data
meas = meas(p,:);
species = species(p);

% gen train and test datasets
half = floor(numObs/2);
training = meas(1:half,:);
trainingSpecies = species(1:half);
sample = meas(half+1:end,:);

% learn on training data and test on sample data
grouphat = classify(sample,training,trainingSpecies);

% Display the confusion matrix for the resulting classification:
group = species(half+1:end);
[C,order] = confusionmat(group,grouphat)




%% Brain data processing
%-- sep_class_112.mat -> brainNet.mat
%-- features_112.mat -> edgeDualNet.mat
%-- names.mat -> cellname stores names of brain cognitive regions
%-- edge_map.mat -> stores mapping between brainNet and edgeDualNet

% load data/edgeDualNet;figure;imshow(abs(S1)*1000,'InitialMagnification',2000);


%-- process brainNet data
clear all;close all; clc;
load data/sep_class_112.mat

data.Xbrain = [];
data.score = [];
data.gnd = [];
data.filename ={};
for i=1:3
    temp = sep_score_class(i).mat;
    nSmp = length(temp);
    Xtmp = zeros(112,112,nSmp);
    Xscore = zeros(1, nSmp);
    filename = cell(1, nSmp);
    for j=1:nSmp
        Xtmp(:,:,j) = temp(j).mat;
        Xscore(j) = temp(j).score;
        filename{j} = temp(j).filename;
    end
    data.Xbrain = cat(3,data.Xbrain,Xtmp);
    data.score = [data.score Xscore];
    
    data.gnd = [data.gnd i*ones(1,nSmp)];
    data.filename = [data.filename filename];
    
end

data.name = 'ADNI brain data';
data.note = 'X edge-dual graphs, Xbrain networks (3 classes 113*133*60';
data.classLabel = {'Normal Control' 'Mild Cognitive Impairment' 'Alzeihmer disease'};


load data/features_112.mat
X = [X1' X2' X3'];

save('data/brain.mat','data');

%-- end process brain data
%%

%-- test converting brainNet to edgedualGraph
clc
n=112;
mat = rand(n,n); % brain net
len = size(mat,1);
nSample = 1;    % no. of brain samples
num_fea = len*(len-1)/2; % no. of edges (from fully connected brainNet)
map = gen_map(len);  % n*n matrix encodes
net = gen_fea_net(map); %-- network of adjacent edges W
X = zeros(nSample, num_fea);
for i=1:nSample
    for k=1:len-1
        for j=k+1:len
            X(i,map(j,k))=mat(j,k);
        end
    end
end

% load 'data/brain.mat'
% data.W = net;
% save('data/brain.mat','data');


% data= genCrossVal(data,5); save ('data/brain.mat','data');
%% END Brain data processing
