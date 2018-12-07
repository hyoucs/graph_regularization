
%% -- DISM with CMUFace straightpose with noise
clear all;close all; clc;
load data/datSDM; 
data = dataAll{5}; 
%%
close all
dataS = data;
%opt.lambda1s = [logspace(2.5, 0, 4) 0]/1000;
opt.lambda1s = [30 20 10];
opt.lambda2s = [20];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = 0.8; opt.verbose = 1;
m = DISMfold(data,opt,1); 
length(find(m.theta~=0))

dataS.gndFea = (m.theta~=0);
opt.gndNet=1; opt.bkImg = 0; ImgNet(dataS,opt);


%% --- PCA-SVM: perform PCA via SVD by keeping .95 

close all; clear all; clc;

load data/syntData.mat
data = dataAll{4}
[U, S, V]=SVD_CutOff(data.X,.95);
data.X = U'*data.X;

mSVM = accuracySVM(data); mSVM.accuCV


%% --- LDA in general case DONE
close all; clear all; clc;

load data/syntData.mat
data = dataAll{4}
opt = []; opt.verbose = 1; %opt.pct = .99
m = LDASVD(data,opt)


%% -- DISM with brain data
clear all;close all; clc;
load data/datSDM; 
data = dataAll{1};
opt.lambda1s = [logspace(2.5, 0, 4) 0]/1000;
opt.lambda2s = [0 10 50 100 200];

opt.bLinear = 0; opt.alpha = .2; opt.k = 50;
opt.nFeaUpd = 1.0; opt.verbose = 1;
tStart = tic; m = DISMfold(data,opt,1); toc(tStart) % m.accuCV

mSVM = accuracySVM(data,4); mSVM.accuCV


fprintf('avgAccuracy: %4.3f (%3.2f)\n',mean(mSVM.accuCV), std(mSVM.accuCV));

    


%%
load data/brain; whos
nClass = length(unique(data.gnd));
for iClass = 3:3%nClass
    load data/brain; whos
%     % remove each class from data
%     clsIdx = find(data.gnd==iClass);
%     data.X(:,clsIdx) =[];
%     data.gnd(:,clsIdx) =[];
%     data.testSet(:,clsIdx) =[];
    
    % groupping classes
    cls1 = (data.gnd == iClass);
    cls2 = ~cls1;
    data.gnd(cls1) = 1;
    data.gnd(cls2) = 2;
    

    data.X = Xnorm(data.X,1);
%     opt.lambda1s = [logspace(2.5, 0, 10) 0]/1000;
%     opt.lambda2s = [0.1 0.5 1 2];
    opt.lambda1s = [0.01];
    opt.lambda2s = [0.5];
    

    opt.bLinear = 0; opt.alpha = .1; opt.k = 50; 
    opt.nFeaUpd = .5; opt.verbose = 1;
    % supLE = supLE(data,opt);
    tStart = tic;m = DISM(data,opt); toc(tStart) % m.accuCV
    mSVM = accuracySVM(data); 
    disp(strcat('avg svm accuracy:',num2str(mean(mSVM.accuCV)), '(' ,num2str(std(mSVM.accuCV)), ')'));
    
end

%%

opt.topFea=find(m.theta~=0);
% subGraphDisp(data.W,opt); 
brainNameDisp(m.theta,data.cellName);




%% ppiData 

clear all;close all; clc;
load 'data/ppiData'; whos

data = dataAll{4};

opt.lambda1s = [logspace(2.5, 0, 4) 0]/1000;
opt.lambda2s = [1, 2]; opt.verbose = 0;
opt.bLinear = 0; opt.alpha = .2; opt.k = 50; 
opt.nFeaUpd = 1; opt.verbose = 1;
opt.type='Eucl'; % remember to update this

m = DISM(data,opt); 





%% -- toyData DONE (note on smp02)
%---------------------------------
%{
clear all;close all; clc;
load 'data/toyData'; whos

gObj = biograph(triu(data.W),data.feaName,'ShowArrows','off') % get(bg2.nodes,'ID')
gObj = view(gObj);

for i=1:size(data.X,1)
    a=corrcoef(data.X(i,:)',data.gnd');
    a(1,2)
end

%% smp01
opt.lambda1s = [0.503];
opt.lambda2s = [0]; opt.verbose = 0;
opt.bLinear = 0; opt.alpha = .35; opt.k = 7; 
opt.nFeaUpd = 1; opt.verbose = 1;
m = DISM(data,opt); m.wAll
figure; stem(abs(m.wAll(:,1)),'Color','k'); xlim([.8 7.2])


%% smp02
opt.lambda1s = [0.667];
opt.lambda2s = [3]; opt.verbose = 0;
opt.bLinear = 0; opt.alpha = .39; opt.k = 7; 
opt.nFeaUpd = 1; opt.verbose = 1;
m = DISM(data,opt); m.wAll
figure; stem(abs(m.wAll(:,1)),'Color','k'); xlim([.8 7.2])

%% smp03
opt.lambda1s = [0.673];
opt.lambda2s = [5]; opt.verbose = 0;
opt.bLinear = 0; opt.alpha = .4; opt.k = 7; 
opt.nFeaUpd = 1; opt.verbose = 1;
m = DISM(data,opt); m.wAll
figure; stem(abs(m.wAll(:,1)),'Color','k'); xlim([.8 7.2])

%% smp04
opt.lambda1s = [0.6];
opt.lambda2s = [2]; opt.verbose = 0;
opt.bLinear = 0; opt.alpha = .4; opt.k = 7; 
opt.nFeaUpd = 1; opt.verbose = 1;
m = DISM(data,opt); m.wAll
figure; stem(abs(m.wAll(:,1)),'Color','k'); xlim([.8 7.2])
%% -- END toyData 
%}


%% -- Gaussian data w/o network


clear all;close all; clc;
load 'data/2Gauss150Dim'; whos
% load 'data/3Gauss150Dim'; whos
% load 'data/dataAllS'; whos
% data = dataAllS{2};



opt.lambda1s = [logspace(2.5, 0, 10) 0]/1000;
% opt.lambda1s = [0.046];
opt.lambda2s = [0];
opt.bLinear = 0; opt.alpha = .1; opt.k = 40; 
opt.nFeaUpd = .5; opt.verbose = 1;
m = DISM(data,opt)



%% ------- multi-class svm test - done
clear all;close all; clc;
load 'data/3Gauss150Dim'; whos
m = accuracySVM(data)
m.accuCV


%% ----------  LDA

clear all;close all; clc;
% load 'data/prostateData'; whos
% load 'data/3Gauss150Dim'; whos
% load 'data/2Gauss150Dim'; whos
% load data/brain; whos
load 'data/dataAllF'; whos
data = dataAll1;
data.X = data.X(1:50,:);
data.X = (data.X - 2*ones(size(data.X)))+rand(size(data.X))/10000;

tic
% data.X = Xnorm(data.X);
[Y,model] = lda(data.X',data.gnd,length(unique(data.gnd))-1);
toc


%%
clear all;close all; clc;
load 'data/toyExamp';
opt.bLinear = 0; opt.alpha = .1; opt.k = 4; 
opt.nFeaUpd = .5; opt.verbose = 1;
m = DISM(data,opt)





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

figure;
subplot(1,2,1); PlotX(data.X,data.gnd,'','Original space',''); grid on;
subplot(1,2,2); PlotX(mappedX,data.gnd,'','Projected space','');



tic
opt.bLinear = 1; opt.alpha = .1; opt.k = 50; opt.lowerBound = 80;
opt.nFeaUpd = .5; opt.verbose = 1;
supLE = supLE(data,opt);
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


% %-- Witten penalizedLDA
% for i=1:length(lambda1s)
%     options.lambda1=lambda1s(i);
%     m = penalizedLDA(data, options, lambda1s(i));
%     w(:,i) = m.v;
% end
% %-- End Witten penalizedLDA


%-- Clemmensen SparseLDA
uniClass = unique(data.gnd);
nClass = length(uniClass);
Y = zeros(nSmp,nClass);
for idxCls=1:nClass
    idxSmp = find(data.gnd==uniClass(idxCls));
    Y(idxSmp,idxCls) = 1;
end
X=data.X';
for i=1:length(lambda1s)
  [beta os] = slda(X,Y, lambda1s(i), 0, 1, 1000);
  w(:,i) = beta;
end
%-- End Clemmensen SparseLDA

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
% %---------------------------- toy example -----------------------------%
% %----------------------------------------------------------------------%
% %-- scenario setup: 7 nodes with 8 subjects
% %-- 1,2,3 are identical to class except last value of node 1
% %-- nodes 6, 7 are identical to 1 yet no edge connecting
% %-- Subnetwork: 1,2,3 as they are strongly connected
% %-- Set W according, strong connection among gndFea, weak connects otherwise 

% data.X = [0 0 1 1 0; 1 0 1 1 0; 0 1 0 0 1; 1 1 0 1 0; 0 0 0 1 1; 1 0 1 1 0; 1 0 1 1 0];
close all; clear all; clc;

data.gnd = [1 0 1 1 0 1 0 0];

data.X = [1 0 1 1 0 1 0 1;
          1 0 1 1 0 1 0 0;
          1 0 1 1 0 1 0 0; 
          0 0 1 0 0 0 1 0;
          0 0 1 0 0 0 1 0;
          1 0 1 1 0 1 0 1;
          1 0 1 1 0 1 0 1];
data.gndFea = [1 1 1 0 0 0 0];      
data.feaName = {'n1', 'n2', 'n3', 'n4', 'n5', 'n6', 'n7'};

data.W = zeros(7,7);
data.W(1,2) = 1;data.W(1,3) = 1; 
data.W(2,4) = 0.5; 
data.W(3,4) = 0.5; 
data.W(4,5) = 1;
data.W(5,6) = 1; data.W(5,7) = 1;

data.W = max(data.W,data.W');

% data.W = [0 1 1 0 0 0 0; 
%           1 0 0 1 0 0 0; 
%           1 0 0 1 0 0 0; 
%           0 1 1 0 1 0 0; 
%           0 0 0 1 0 1 1; 
%           0 0 0 0 1 0 0; 
%           0 0 0 0 1 0 0];

gObj = biograph(triu(data.W),data.feaName,'ShowArrows','off');
gtNodeIdx= find(data.gndFea==1); 
set(gObj.nodes(gtNodeIdx),'Color',[1 0 0],'size',[40 30]);
gObj = view(gObj);
save('data/toyData.mat','data');




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


%% ------- Simple Network Visualization
close all; clear all; clc
load data/simpNet

gObj = biograph(triu(data.W),data.feaName,'ShowArrows','off') % get(bg2.nodes,'ID')
gObj = view(gObj);

%-- color GTnodes
gtNodeIdx= find(data.gndFea==1); %-- find gtNodes in top 100 nodes
set(gObj.nodes(gtNodeIdx),'Color',[1 0 0],'size',[40 30]);
dolayout(gObj);







%% -- Plot lines connecting points

close all;
x = [0 1 1 0 0 1; ...
     1 1 0 0 0 0];
y = [0 0 1 1  0 0; ...
     0 1 1 0 1 1 ];
 
h = figure;
for i=1:size(x,2)
    line([x(1,i) y(1,i)],[x(2,i) y(2,i)],'Marker','.','LineStyle','-');
    figure(h);
end
xlim([-.2 1.2]);ylim([-.2 1.2]);



%-------------------------------------------------------------------------%
%------------------------------ CMU preprocessing-------------------------%
%-------------------------------------------------------------------------%
%{

%% -- step 1: put same-size images with same size from each folder (forget
% this step after having data in mat format)

% close all; clear all; clc;
% in_folder = 'data/faces';
% 
% 
% folds = dir(strcat(in_folder,''));
% folds(1:3,:) = []; % remove 3 first folders
% 
% out_folder = 'data/ALL60x64';
% for j=1:size(folds,1)
%     stm=['mv ',in_folder,'/', folds(j).name,'/*2.pgm ',out_folder];
%     dos(stm);
% end
% disp('done');
% 
% out_folder = 'data/ALL30x32';
% for j=1:size(folds,1)
%     stm=['mv ',in_folder,'/', folds(j).name,'/*4.pgm ',out_folder];
%     dos(stm);
% end
% disp('done');
% 
% out_folder = 'data/ALL120x128';
% for j=1:size(folds,1)
%     stm=['mv ',in_folder,'/', folds(j).name,'/*.pgm ',out_folder];
%     dos(stm);
% end
% disp('done');
%}

%{
%% %-- step 2: convert into mat format (after this, forget step 1 forever!)
% 
% close all; clear all; clc;
% % in_folder = 'data/ALL30x32';
% % in_folder = 'data/ALL60x64';
% in_folder = 'data/ALL120x128';
% 
% files = dir(strcat(in_folder,''));
% files(1:2,:) = []; % remove 2 first folders
% nFile = size(files,1);
% 
% I = imread([in_folder '/' files(1).name]);
% [r,c] = size(I);
% 
% X = zeros(r*c,nFile);
% 
% 
% gndPerson = zeros(1, nFile);
% gndPosition = zeros(1, nFile);
% gndExpression = zeros(1, nFile);
% gndGlass = zeros(1, nFile);
% 
% for j=1:nFile
%     I = imread([in_folder '/' files(j).name]);
%     X(:,j)= I(:);
%     if (findstr(files(j).name,'an2i')>0)     c1 = 1; end;
%     if (findstr(files(j).name,'at33')>0)     c1 = 2; end;
%     if (findstr(files(j).name,'boland')>0)   c1 = 3; end;
%     if (findstr(files(j).name,'bpm')>0)      c1 = 4; end;
%     if (findstr(files(j).name,'ch4f')>0)     c1 = 5; end;
%     if (findstr(files(j).name,'cheyer')>0)   c1 = 6; end;
%     if (findstr(files(j).name,'choon')>0)    c1 = 7; end;
%     if (findstr(files(j).name,'danieln')>0)  c1 = 8; end;
%     if (findstr(files(j).name,'glickman')>0) c1 = 9; end;
%     if (findstr(files(j).name,'karyadi')>0)  c1 = 10; end;
%     if (findstr(files(j).name,'kawamura')>0) c1 = 11; end;
%     if (findstr(files(j).name,'kk49')>0)     c1 = 12; end;
%     if (findstr(files(j).name,'megak')>0)    c1 = 13; end;
%     if (findstr(files(j).name,'mitchell')>0) c1 = 14; end;
%     if (findstr(files(j).name,'night')>0)    c1 = 15; end;
%     if (findstr(files(j).name,'phoebe')>0)   c1 = 16; end;
%     if (findstr(files(j).name,'saavik')>0)   c1 = 17; end;
%     if (findstr(files(j).name,'steffi')>0)   c1 = 18; end;
%     if (findstr(files(j).name,'sz24')>0)     c1 = 19; end;
%     if (findstr(files(j).name,'tammo')>0)    c1 = 20; end;
%     
%     if (findstr(files(j).name,'straight')>0) c2 = 1; end;
%     if (findstr(files(j).name,'left')>0)     c2 = 2; end;
%     if (findstr(files(j).name,'right')>0)    c2 = 3; end;
%     if (findstr(files(j).name,'_up_')>0)     c2 = 4; end;
%     
%     if (findstr(files(j).name,'neutral')>0)  c3 = 1; end;
%     if (findstr(files(j).name,'happy')>0)    c3 = 2; end;
%     if (findstr(files(j).name,'sad')>0)      c3 = 3; end;
%     if (findstr(files(j).name,'angry')>0)    c3 = 4; end;
% 
%     if (findstr(files(j).name,'open')>0)        c4 = 1; end;
%     if (findstr(files(j).name,'sunglasse')>0)   c4 = 2; end;
% 
%     gndPerson(j)=c1;
%     gndPosition(j)=c2;
%     gndExpression(j)=c3;
%     gndGlass(j)=c4;
% 
% end
% 
% % X30 = X;
% % X60 = X;
% X120 = X;
% 
% data.X30 = X30;
% data.X60 = X60;
% data.X120 = X120;
% data.gndPerson = gndPerson;
% data.gndPosition = gndPosition;
% data.gndExpression = gndExpression ;
% data.gndGlass = gndGlass;
% 
% data.perClTxt='1-20: an2i, at33, boland, bpm, ch4f, cheyer, choon, danieln, glickman, karyadi, kawamura, kk49, megak, mitchell, night, phoebe, saavik, steffi, sz24, and tammo';
% data.posClTxt='1-4: straight, left, right, up';
% data.expClTxt='1-4: neutral, happy, sad, angry';
% data.glaClTxt='1-2: open, sunglasses';
% data.rc30 = [30 32];
% data.rc60 = [60 64];
% data.rc120 = [120 128];
% 
% save('data/CMUallResolutions.mat','data')
%}

%{


%% -- step 3: Test
close all; clear all; clc;
load 'data/CMUallResolutions.mat'

i=10; 
r = data.rc30(1); c = data.rc30(2);
subplot(1,3,1); imshow(reshape(data.X30(:,1),r,c),[]); 
r = data.rc60(1); c = data.rc60(2);
subplot(1,3,2); imshow(reshape(data.X60(:,1),r,c),[]); 
r = data.rc120(1); c = data.rc120(2);
subplot(1,3,3); imshow(reshape(data.X120(:,1),r,c),[]); 

%}

%{


%% -- step 4: Extract 2 classes: straight pose w and w/o glasses

close all; clear all; clc;
load 'data/CMUallResolutions.mat'
dataFull = data; clear data;

% data.X = dataFull.X60;            %-- all poses
% data.gnd = dataFull.gndGlass;


% idx = find(dataFull.gndPosition == 1);   %-- straigth pose images
idx = find(dataFull.gndPosition ~= 4);   %-- all poses except up one
data.X = dataFull.X60(:,idx);            %-- get from 60x64 resolution
data.gnd = dataFull.gndGlass(:,idx);


r = dataFull.rc60(1); 
c = dataFull.rc60(2); 

%--visualize selected images
nImg=size(data.X,2);
nRow=10;
figure; fpr=ceil(nImg/nRow);
for i=1:nImg
    subplot(nRow,fpr,i); imshow(reshape(data.X(:,i),r,c),[]); 
end

data.r = r;
data.c = c;

%--visualize wearing glasses images
idx=find(data.gnd == 2);  
r = data.r; c = data.c; 
Xglass = data.X(:,idx);
figure; nImg=size(Xglass,2);nRow=10;
fpr=ceil(nImg/nRow);
for i=1:nImg
    subplot(nRow,fpr,i); imshow(reshape(Xglass(:,i),r,c),[]); 
end


muX1 = mean(data.X(:,find(data.gnd == 1))')';
muX2 = mean(data.X(:,find(data.gnd == 2))')';
figure;
subplot(1,2,1);imshow(reshape(muX1,r,c),[]); 
subplot(1,2,2);imshow(reshape(muX2,r,c),[]); 
data.clsMean = [muX1 muX2];

% save('data/CMU60GlassAllPose.mat','data');

save('data/CMU60Glass.mat','data');
%}

%{


%% -- step 5: building network and downsample to nodes

close all; clear all; clc;
load 'data/CMU60Glass.mat';

pctConst = .5; %-- percentage of pixels kept from an image
k = 5; %-- connect to k neighbors from each nodes

%-- generate data points based on locations 
syntXfull = zeros(2,data.r*data.c);
idx = 0;
for i=1:data.c
    for j=-1:-1:-data.r %-- resemble an image indexed from TOP-LEFT corner
%     for j=1:data.r
        idx = idx+1;
        syntXfull(:,idx) = [i,j]';
    end
end

%-- keep pctConst points
idxSet = randperm(length(syntXfull));
idxSele = idxSet(1:ceil(pctConst*length(syntXfull)));
syntX = syntXfull(:,idxSele);


%-- Compute Affinity matrix for network topo
opt = [];
opt.type = 'Eucl';
[W] = PWdistance(syntX,syntX,opt);
[dmp idx] = sort(W,2,'ascend'); 

idx(:,[1,k+2:end]) = []; %-- remove 1st (selfnode sim) + last smallest sim
kNN = zeros(size(W));
nSmp = size(syntX,2);
for i = 1 : nSmp
    kNN(i,idx(i,:)) = 1;
end
W = W.*kNN;

%-- normalize W: set smallest non-zero = 1 due to Eucl distance
W = W/max(max(W));
B = W; B(~B) = inf;
W = 1 + min(min(B)) - W;
W = W.*kNN;
W = max(W,W');  

% -- visualize network on an image sample.
halfW = triu(W);
[row col] = find(halfW>0);

x = syntX(:,row);
y = syntX(:,col);

figure;
h1 = subplot(1,3,1); PlotX(syntX,'','','','');
h2 = subplot(1,3,2); imshow(reshape(data.clsMean(:,2),data.r,data.c),[]); hold on
imgX = ones(data.r*data.c,1)*256;
imgX(idxSele) = data.clsMean(idxSele,2);
h3 = subplot(1,3,3);imshow(reshape(imgX,data.r,data.c),[]); 

for i=1:size(x,2)
    subplot(h1);
    line([x(1,i) y(1,i)],[x(2,i) y(2,i)],'LineStyle','-');
    subplot(h2);
    line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'Marker','.','LineStyle','-');
%     subplot(h3);
%     line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'Marker','.','LineStyle','-');
end


data.Xfull = data.X;
data.X = data.X(idxSele,:);
data.Wnorm = W;
data.W = W>0;
data.syntXfull = syntXfull;
data.idxSele = idxSele;
data.syntX = syntX;

data.name = 'CMUFace straight pose w and w/o glasses';
data = genCrossVal(data,5);

save('data/CMU60Glass.mat','data');

%}

%{



%% -- step 6: check
close all; clear all; clc;
% load 'data/CMU60Glass.mat';
% 0->black, 255-> white
load data/datSDM.mat; 
data = dataAll{3}; % front pose
 

halfW = triu(data.W);
[row col] = find(halfW>0);

x = data.syntX(:,row);
y = data.syntX(:,col);

img = data.clsMean(:,2);
% img = data.Xfull(:,4);


figure;
h1 = subplot(1,3,1); PlotX(data.syntX,'','','','');
h2 = subplot(1,3,2); imshow(reshape(img,data.r,data.c),[]); hold on

imgX = ones(data.r*data.c,1)*255;
imgX(data.idxSele) = img(data.idxSele);
% imgX(data.idxSele) = data.clsMean(data.idxSele,2);
h3 = subplot(1,3,3);imshow(reshape(imgX,data.r,data.c),[]); 

for i=1:size(x,2)
    subplot(h1);
    line([x(1,i) y(1,i)],[x(2,i) y(2,i)],'LineStyle','-');
    subplot(h2);
    line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    subplot(h3);
    line([x(1,i) y(1,i)],[x(2,i)*(-1) y(2,i)*(-1)],'LineStyle','-');
    
end

%}

%{



%% -- step 7: Manually generate gndFea for CMUface
close all; clear all; clc;
% 0->black, 255-> white
load data/datSDM.mat; 
data = dataAll{3}; % front pose

halfW = triu(data.W);
[row col] = find(halfW>0);

x = data.syntX(:,row);
y = data.syntX(:,col);

img = data.clsMean(:,2);
% img = data.Xfull(:,4);


close all
%----------------- BEGIN find gndFea for front pose  ---------------------%
img = reshape(img,data.r,data.c);
rrng = (24:31); crng = (20:43); % front pose

% % create a rectangular in the real img:
% img(rrng(1),crng) = 255; img(rrng(end),crng) = 255;
% img(rrng,crng(1)) = 255; img(rrng,crng(end)) = 255;

% create a corresponding rectangular block:
bimg = 0*ones(data.r,data.c);
bimg(rrng,crng) =1;
bimg = bimg.*img;

img = img(:);
bimg = bimg(:);
bimg = img.*(bimg<1.4*mean(img) & bimg >0);
dataAll{3}.gndFea = bimg(data.idxSele)>0; % front pose
% dataAll{2}.gndFea = bimg(data.idxSele)>0; % all poses
length(find(dataAll{3}.gndFea)) % 79 = number of gndFea



figure;
% h1 = subplot(1,3,1); PlotX(data.syntX,'','','','');
h1 = subplot(1,3,1); imshow(reshape(bimg,data.r,data.c),[]); 
h2 = subplot(1,3,2); imshow(reshape(img,data.r,data.c),[]); hold on

imgX = ones(data.r*data.c,1)*255;
imgX(data.idxSele) = img(data.idxSele);
% imgX(data.idxSele) = data.clsMean(data.idxSele,2);
h3 = subplot(1,3,3);imshow(reshape(imgX,data.r,data.c),[]); 
%----------------- END find gndFea for front pose  -----------------------%



%%
close all
%----------------- BEGIN find gndFea for ALL poses  ---------------------%
data = dataAll{2}; % all poses
% gndPose = 1; %-- 1-straight 2-right 3-left 4-up
% idxPose = find(data.gndPosition==gndPose);
idxGlass = find(data.gnd == 2);
% idxGlass = idxGlass(ismember(idxGlass,idxPose));

img = mean(data.Xfull(:,idxGlass),2);
% img = data.Xfull(:,4);

img = reshape(img,data.r,data.c);
rrng1 = (24:31); crng1 = (15:47); %-- all except up-pose
rrng2 = (16:22); crng2 = (19:43); %-- up-pose

% % create a rectangular in the real img:
% img(rrng1(1),crng1) = 255; img(rrng1(end),crng1) = 255;
% img(rrng1,crng1(1)) = 255; img(rrng1,crng1(end)) = 255;

% create a corresponding rectangular block:
bimg = 0*ones(data.r,data.c);
bimg(rrng1,crng1) =1;
bimg(rrng2,crng2) =1;
bimg = bimg.*img;

img = img(:);
bimg = bimg(:);
bimg = img.*(bimg<1.6*mean(img) & bimg >0);
dataAll{2}.gndFea = bimg(data.idxSele)>0; % all poses


figure;
% h1 = subplot(1,3,1); PlotX(data.syntX,'','','','');
h1 = subplot(1,3,1); imshow(reshape(bimg,data.r,data.c),[]); 
h2 = subplot(1,3,2); imshow(reshape(img,data.r,data.c),[]); hold on

imgX = ones(data.r*data.c,1)*255;
imgX(data.idxSele) = img(data.idxSele);
% imgX(data.idxSele) = data.clsMean(data.idxSele,2);
h3 = subplot(1,3,3);imshow(reshape(imgX,data.r,data.c),[]); 
%----------------- END find gndFea for ALL poses  -----------------------%

%}

%{

%% -- step 8: Visualize subnetwork 
close all; clear all; clc;
load data/datSDM.mat; 
data = dataAll{3}; % front pose
opt.gndNet=1; opt.bkImg = 0; ImgNet(data,opt);

%----------------------- END CMU preprocessing --------------------------%
%}

%-------------------------------------------------------------------------%
%% -- Prepare Brain (delete MCI class) and CMU face data

%{

close all; clear all; clc;
load data/brain.mat
hist1D(data.gnd);

%-- del MCI 2nd class
idx = find(data.gnd==2);
data.X(:,idx) = [];
data.gnd(idx) = [];
data.score(idx) = [];
data.filename(idx) = [];
data.name = 'ADNI brain data 2 classes NC n AD';
data.classLabel(2) = [];
data.note = 'X: edge-dual graphs; Xbrain networks 2 classes 113/60;'
data.Xbrain(:,:,idx) = [];
data.testSet(idx) = [];

dataAll = {};
dataAll{1} = data;


%-- del AD 3rd class
idx = find(data.gnd==3);
data.X(:,idx) = [];
data.gnd(idx) = [];
data.score(idx) = [];
data.filename(idx) = [];
data.name = 'ADNI brain 2 classes NC n MCI';
data.classLabel(2) = [];
data.note = 'X: edge-dual graphs; Xbrain networks 2 classes 113/133;'
data.Xbrain(:,:,idx) = [];
data.testSet(idx) = [];


load data/datSDM
dataAll{end+1} = data;
save('data\datSDM','dataAll');


load 'data/CMU60GlassPos123.mat'
dataAll{2} = data;
load 'data/CMU60Glass.mat';
dataAll{3} = data;
save('data/datSDM.mat','dataAll');

%% -- End prepare Brain and CMU face data
%}

