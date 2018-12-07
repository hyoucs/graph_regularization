function [m] = DISM(data,opt)
%------------------------------------------------------------------------%
% DIscriminant Subnetwork Feature Mining from supervised network data
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
%           .lambda2s:   set of L2 tradeoff for smoothness
%           .minFea: smallest number of selective features for each embedded dim
%           .verbose: print results on-the-fly
% Output:
%     + model: struct data type
%           .Y: X in embedding space
%           .Z: X in transformed space
%           .theta: coeff to form each dim of transformed space
%           .predictedLabel: labels predicted by SVM on all nFold testdata
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

data.X = Xnorm(data.X);
[nFea nSmp] = size(data.X);


nFold = 1;
if isfield(data,'testSet'),
    nFold = length(unique(data.testSet));
end

verbose=0;
if isfield(opt,'verbose'),
    verbose= opt.verbose;
end


if isfield(opt,'lambda1s'),
    lambda1s = opt.lambda1s;
else
    lambda1s = [logspace(2.5, 0, 4) 0];
    lambda1s = lambda1s/1000;
end


if isfield(opt,'lambda2s'),
    lambda2s = opt.lambda2s;
else
    lambda2s = sort([logspace(.5, 0, 3) .5 .01],'ascend');
end

%(1) Find X in embedding space
if isfield(data,'Y')
    Y = data.Y;
else
    % not provided, find Y by supLE
    sModel = supLE(data,opt);
    Y = sModel.Y;
    m.Y = sModel.Y;
    m.eVecs = sModel.eVecs;
    m.eVals = sModel.eVals;
end

[dDim nSmp] = size(Y);

%(2) regression embedding data Y

if verbose
    fprintf('%3s %6s %6s %11s %5s %5s %6s \n','Idx','lmbd1','lmbd2','cntFea','avgAC','stdAC','elpTime');
end

% normalized features, test without this later!
for iFea = 1:dDim
    Y(iFea,:) = Xnorm(Y(iFea,:));
end

%     % wLS found from least square regression
%     wLS = data.X'\data.y';
%     denom = sum(abs(wLS'));

% store sparse coeff vectors into w
m.wAll = zeros(nFea, dDim*length(lambda1s)*length(lambda2s));
w = zeros(nFea, dDim);

% accuCV store L1, L2 and nFold accracies
m.accuCV = zeros(length(lambda1s)*length(lambda2s),2+nFold);

dataProj.gnd = data.gnd;
if isfield(data,'testSet')
    dataProj.testSet = data.testSet;
end
bestIdx = 1;
bestAcc = 0.1;
idx = 0;
for i = 1:length(lambda1s)
    opt.lambda1 = lambda1s(i);
    for j = 1:length(lambda2s)
        opt.lambda2 = lambda2s(j);
        idx = idx +1;

        tic;
        % regress each dim of embedded data
        for iFea = 1:dDim
            data.y = Y(iFea,:);
            regModel = regCoordAscent(data, opt);
            w(:,iFea) = regModel.a;
        end
        
        cntFea = sum((w(:,:)~=0));
        
        if isfield(data,'testSet')
            
          
            dataProj.X = w'*data.X;
            svmModel = accuracySVM(dataProj);
            m.accuCV(idx,1:2) = [opt.lambda1 opt.lambda2]; 
            m.accuCV(idx,3:end) = svmModel.accuCV;
            %m.predictedLabel = svmModel.predictedLabel(:,1);
            
            % find best lambda1, lambda2 given limited features
            currAcc =  mean(svmModel.accuCV);
            if (currAcc>bestAcc) & (cntFea<=100)
                bestIdx = idx;
                bestAcc = currAcc;
            end
        end
    
        if verbose
            fprintf('%3d %6.3f %6.3f %11s %5.2f %5.2f %6.2f \n',i,lambda1s(i),lambda2s(j),num2str(cntFea),mean(m.accuCV(idx,3:end)),std(m.accuCV(idx,3:end)),toc);
        end
        m.wAll(:,dDim*(idx-1)+1:dDim*idx ) = w;
       
    end    
end


% save variables to output model
% m.nFeaUpd = regModel.nFeaUpd;
m.lambda1 = m.accuCV(bestIdx,1);
m.lambda2 = m.accuCV(bestIdx,2);
m.acc = bestAcc;
m.theta = m.wAll(:,dDim*(bestIdx-1)+1:dDim*bestIdx);
m.Z = m.theta'*data.X;


% % plot proj space and selected features
% figure;
% nModels = size(m.accuCV,1);
% for i = 1:nModels
%     theta = m.wAll(:,dDim*(i-1)+1:dDim*i);
%     Z = theta'*data.X;
%     subplot(2,nModels,i);PlotX(Z,data.gnd,'','projSpace',''); 
%     subplot(2,nModels,i + nModels);stem(abs(theta(:,1)),'Color','k'); title(strcat('selFea',num2str(mean(m.accuCV(i,3:end)),2))); xlim([.8 7.2]);
%     if (dDim > 1) % plot only the 2nd one
%         hold on;
%         stem(abs(theta(:,2)),'Color','r')
%     end
% end



% %-- merging all dDim in theta here:
% % figure;PlotX(m.Y,data.gnd,'','embedded space',''); 
% %     % visualize pathways of coefficients
% %     s = sum(abs(w'),2)/denom;
% %     figure;plot(s, w', '-o'); title(strcat('lambda2 = ',num2str(opt.lambda2)));



end