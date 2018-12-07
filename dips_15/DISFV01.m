function [m] = DISF(data,opt)
% function [m] = feaSubnetV01(data,opt)
%------------------------------------------------------------------------%
% Find optimal subnets from supervised network data
%  accuracy computed separately from lambda1
% 2 steps: 
%       (1) (nonlinear)supervisedLE; 
%           skip (1) if data.Y is provided 
%       (2) (linear)regCoordAscent
% 
% Input:
%     + data: struct data type
%           .X [nFea nSmp]: dataset X (no need?)
%           .gnd [1 nSmp]: class labels
%           .W [nFea nFea]: network topo
%           .lambda1s [1 n]: set of lambda1
%           .Y [mFea nSmp]: X in embedding space, if not provided, compute
%           via supervisedLE
%     + opt: structure
%           .bLinear: linear or non-linear LE (default 0-nonlinear)
%           .alpha: tradeoff btw Lb and Aw
%           .k:   number of nearest neighbors
%           .lambda1s:  set of L1 tradeoff for sparsness 
%           .lambda2:   L2 tradeoff for smoothness
%           .lowerBound: number of selective features for each embedded dim
%           .verbose: print on-the-fly steps
% Output:
%     + model: struct data type
%           .Y: X in embedding space
%           .Z: X in transformed space
%           .beta: coeff to form each dim of transformed space
%           .predictedLabel: labels predicted by SVM on all nFold testdata
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%



data.X = Xnorm(data.X);
[nFea nSmp] = size(data.X);

lowerBound = 10;
if isfield(opt,'lowerBound'),
    lowerBound = opt.lowerBound;
end


if ~isfield(opt,'lambda2'),
    opt.lambda2 = 1;
end


verbose=0;
if isfield(opt,'verbose'),
    verbose= opt.verbose;
end

if isfield(opt,'lambda1s'),
    lambda1s = opt.lambda1s;
else
    lambda1s = [logspace(2.5, 0, 40) 0]; 
    lambda1s = lambda1s/100;
end



%(1) Find X in embedding space
if isfield(data,'Y') 
    Y = data.Y;
else
    % not provided, find Y by supervisiedLE
    supLE = supervisedLE(data,opt);
    Y = supLE.Y;
    m.Y = supLE.Y;
    m.eVecs = supLE.eVecs;
    m.eVals = supLE.eVals;
end


[dDim nSmp] = size(Y);

%(2) regression embedding data

if verbose
    fprintf('%7s %10s %8s %10s \n','lmbdIdx','lmbd1','cntFea','elpTime');
end

beta = zeros(nFea, dDim);
for iFea = 1:dDim % regress each dim of embedded data
    data.y = Xnorm(Y(iFea,:));
    % beta wLS found from least square regression
    wLS = data.X'\data.y';
    denom = sum(abs(wLS'));
    
    w = zeros(nFea, length(lambda1s));
    for i = 1:length(lambda1s)
        tic;
        opt.lambda1 = lambda1s(i);
        regModel = regCoordAscent(data, opt); %-- my own
        w(:,i) = regModel.a;
        %     w(:,i) = eNet(data,options); %-- elasticNet with network
        cntFea = sum((w(:,i)~=0));
        if verbose
            fprintf('%7d %10.4f %8d %10.2f \n',i,lambda1s(i),cntFea,toc);
        end
    end
    if verbose, fprintf('\n');end

    
    % visualize pathways of coefficients
    s = sum(abs(w'),2)/denom;
    figure;plot(s, w', '-o'); title(strcat('lambda2 = ',num2str(opt.lambda2)));
    
    % choose best w to save to beta
    bestIdx = 1;
    for wIdx = 2:size(w,2)
        cntNotZero = find(w(:,wIdx)~=0);
        if length(cntNotZero) > lowerBound
            bestIdx=wIdx;
            break;
        end
    end
    
    beta(:,iFea) = w(:,bestIdx);
end

% save variables to output model
m.nFeaUpd = regModel.nFeaUpd;
m.lambda1 = lambda1s(bestIdx);
m.lambda2 = opt.lambda2;
m.lowerBound = lowerBound;
m.beta = beta;
m.Z = beta'*data.X;

figure; 
subplot(2,2,1);PlotX(m.Y,data.gnd,'','',''); title('embbed space');
subplot(2,2,2);PlotX(m.Z,data.gnd,'','',''); title('projected space');

subplot(2,2,3);stem(abs(m.beta(:,1))); title('1st principle');
if size(beta,2) > 1
    subplot(2,2,4);stem(abs(m.beta(:,2))); title('2nd principle');
end

%(3) evaluate ACCURACY via nFold-CV

if isfield(data,'testSet') 
    nFold = length(unique(data.testSet));
    m.accuCV = zeros(1,nFold);
    m.outLabel = zeros(1,nSmp);
    
    for iFold = 1:nFold 
        % setup test and train data
        testIdx = (data.testSet == iFold);
        trainIdx = ~testIdx;
        
        trainData.X = m.Z(:,trainIdx);
        trainData.gnd = data.gnd(trainIdx);
        
        testData.X = m.Z(:,testIdx);
        testData.gnd = data.gnd(testIdx);
       
        outlabel = multisvm(trainData.X',trainData.gnd',testData.X');
        m.accuCV(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
        m.predictedLabel(find(data.testSet == iFold)) = outlabel;
    end 
    m.predictedLabel = [data.gnd' m.predictedLabel'];
end

end