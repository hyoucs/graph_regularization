function [m, m2] = dipsFold(data,opt,idxFold)
%------------------------------------------------------------------------
% Learning discriminant subnetwork from supervised network data
% 2 steps:
%       (1) (nonlinear) supervisedLE;
%           skip (1) if data.Y is provided
%       (2) (linear) regCoordAscent
%
% Input:
%     + data: struct data type
%           .X [nSmp nFea]:     dataset X
%           .gnd [nSmp 1]:      class labels
%           .W [nFea nFea]:     network topo
%           .lambda1s [1 n]:    set of lambda1
%           .Y [nSmp mFea]:     X in embedding space
%                               (if not provided, compute via supervisedLE)
%           .A [nSmp nSmp]:     pairwise similarity (among samples)
%     + opt: structure
%           .bLinear:           linear or non-linear LE (default 0-nonlinear)
%           .beta:             tradeoff btw Lb and Aw
%           .k:                 number of nearest neighbors
%           .lambda1s:          set of L1 tradeoff for sparsness
%           .lambda2s:          set of L2 tradeoff for smoothness
%           .minFea:            smallest number of selective features for each embedded dim
%           .verbose:           print results on-the-fly
%     + idxFold: compute accuracy on this fold only
% Output:
%     + model: struct data type
%           .Y:                 X in embedding space
%           .Z:                 X in transformed space
%           .theta:             coeff to form each dim of transformed space
%           .predictedLabel:    labels predicted by SVM on all nFold testdata
%------------------------------------------------------------------------







    % - - - - - - - - - - - - - IMPORT PARAMETERS - - - - - - - - - - - - - - - - - 

    verbose=0;
    if isfield(opt,'verbose')
        verbose= opt.verbose;
    end

    if isfield(opt,'lambda1s')
        lambda1s = opt.lambda1s;
    else
        lambda1s = [logspace(2.5, 0, 4) 0];
        lambda1s = lambda1s/1000;
    end

    if isfield(opt,'lambda2s')
        lambda2s = opt.lambda2s;
    else
        lambda2s = sort([logspace(.5, 0, 3) .5 .01],'ascend');
    end

    % For the final step classification
    % if svm is not explicitly required, set as required by default
    if ~isfield(opt, 'svm')
        opt.svm = true;
    end






    % - - - - - - - - - - - - - IMPORT CV TRAINING SETTINGS - - - - - - - - - - - - - - - - - 


    % - - - - - - Determine # (cv folds)
    % data.testSet: indices of test samples for each cv fold
    if ~isfield(data,'testSet')
        error('This function requires cross validation indexing!');
    else
        nFold = length(unique(data.testSet));
    end

    % - - - - - - Prepare training samples and labels
    % normalized all features
    data.X = Xnorm(data.X')';
    % if idxFold is not provided,
    trainData = data;
    % if idxFold (3rd parameter of this function) is provided, 
    % then perform training and testing on a single fold indicated by idxFold.
    trainIdx    = (trainData.testSet ~= idxFold);
    trainData.X = trainData.X(trainIdx,:);
    trainData.gnd = trainData.gnd(trainIdx,1);
    testIdx     = (data.testSet == idxFold);
    testData.gnd = data.gnd(testIdx,1); 
    if isfield(trainData,'Y')
        trainData.Y = data.Y(trainIdx,:);
    end
    % trainData.X = Xnorm(trainData.X')';
    [nSmp, nFea] = size(trainData.X);

    % - - - - - - Prepare truncated pairwise similarity matrix
    trainData.A = data.A(trainIdx, trainIdx);







    % - - - - - - - - - - - - - DO WORK - - - - - - - - - - - - - - - - - 

    % - - - - - - STEP 1: Find X in embedding space - - - - - - - - - - -
    if isfield(trainData,'Y') 
        % embbedding space Y is provided
        Y = trainData.Y;
    else    
        % embedding space Y is not provided, find it by dipsLE.m
        sModel = dipsLE(trainData, opt);
        Y = sModel.Y;
        % save embedding space information to structure
        m.Y = Y;
        m.eVecs = sModel.eVecs;
        m.eVals = sModel.eVals;
    end
    % visualProj(Y, trainData.gnd, ['trainData Fold #', num2str(idxFold)]);
    [nSmp, d] = size(Y);

    % if svm is suppressed, then will use kmeans instead of linear svm
    % the center of each class should be computed once you get Y in the
    % embedding space (d-dimension) limited to training samples. 
    if ~opt.svm
        m.gnd = trainData.gnd;
        m = classCenter(m);
    end
    
    % % plot subjects in embedding space;
    % if opt.bLinear 
    %     Xproj = data.X*Y;
    %     figure; subplot(1,2,1); gscatter(Xproj(trainIdx,1),Xproj(trainIdx,2),trainData.gnd);
    %     subplot(1,2,2); gscatter(Xproj(testIdx,1),Xproj(testIdx,2),testData.gnd);
    % else
    %     figure; gscatter(Y(:,1),Y(:,2),trainData.gnd);
    % end

    % figure;
    % class_idx = find(trainData.gnd == 1);
    % plot3(Y(class_idx,1),Y(class_idx,2),Y(class_idx,3), 'r*');
    % hold on;  figure;
    % gscatter(Y(:,2),Y(:,3),trainData.gnd);
    % class_idx = find(trainData.gnd == 2);
    % plot3(Y(class_idx,1),Y(class_idx,2),Y(class_idx,3), 'b*');
    % grid on;
    

    % - - - - - - STEP 2: regression Y - - - - - - - - - - - - - - - - -
    % table column header
    if verbose
        fprintf('%3s %6s %6s %11s %5s %5s %6s \n',...
            'Idx','lmbd1','lmbd2','cntFea','avgAC','stdAC','elpTime');
    end

    % normalized features, test without this later!
    for iFea = 1:d
        Y(:,iFea) = Xnorm(Y(:,iFea)')';
    end

    % store parameter vectors into w
    m.UAll = zeros(nFea, d*length(lambda1s)*length(lambda2s));
    U = zeros(nFea, d);

    % dataProj: dataset in embedding space
    dataProj.gnd = data.gnd;
    dataProj.testSet = data.testSet;

    if ~opt.svm
        dataProj.classCenter = m.classCenter;
    end

    % record best parameter setting corresponding to highest accuracy 
    bestIdx = 1;
    bestAcc = 0.1;
    idx = 0;    % - - index for parameter combination
    m2 = {};
    for i = 1:length(lambda1s) % - - loop all sparse parameters
        opt.lambda1 = lambda1s(i);
        
        for j = 1:length(lambda2s)  % - - loop all network constraint parameters
            opt.lambda2 = lambda2s(j);
            idx = idx + 1;   

            % regress each dim of embedded data Y,
            % call dipsCoordAscent
            tic;
            for iFea = 1:d
                trainData.y = Y(:, iFea);
                regModel = dipsCoordAscent(trainData, opt);
                U(:,iFea) = regModel.u;
            end

            % count non-zero features for all dimensions
            cntFea = sum((U(:,:)~=0));

            % project dataset with learned sparse coefficients
            dataProj.X = data.X*U; % proj entire dataset
            % visualProj(dataProj.X, dataProj.gnd, ['testData Fold #', num2str(idxFold)]);

            % classify test samples with learned model
            svmModel = accuracySVM(dataProj, idxFold);
            kmsModel = accuracyKmeans(dataProj, idxFold);

            % export models 
            m2{i}{j}.lambda1 = lambda1s(i);
            m2{i}{j}.lambda2 = lambda2s(j);
            m2{i}{j}.accuCV(idx,1:2) = [opt.lambda1 opt.lambda2];
            m2{i}{j}.accuCV(idx,3) = svmModel.accuCV;
            m2{i}{j}.accuCV(idx,4) = kmsModel.accuCV;
            m2{i}{j}.UAll(:,d*(idx-1)+1:d*idx ) = U;
            %m.predictedLabel = svmModel.predictedLabel(:,1);

            if verbose
                fprintf('%3d %6.3f %6.3f %11s %5.2f %5.2f %6.2f \n',...
                    i,lambda1s(i),lambda2s(j),num2str(cntFea),...
                    svmModel.accuCV, kmsModel.accuCV, toc);
            end

        end

    end


    % % plot proj space and selected features
    % figure;
    % nModels = size(m.accuCV,1);
    % for i = 1:nModels
    %     theta = m.wAll(:,d*(i-1)+1:d*i);
    %     Z = theta'*data.X;
    %     subplot(2,nModels,i);PlotX(Z,data.gnd,'','projSpace','');
    %     subplot(2,nModels,i + nModels);stem(abs(theta(:,1)),'Color','k'); title(strcat('selFea',num2str(mean(m.accuCV(i,3:end)),2))); xlim([.8 7.2]);
    %     if (d > 1) % plot only the 2nd one
    %         hold on;
    %         stem(abs(theta(:,2)),'Color','r')
    %     end
    % end



    % %-- merging all d in theta here:
    % % figure;PlotX(m.Y,data.gnd,'','embedded space','');
    % %     % visualize pathways of coefficients
    % %     s = sum(abs(w'),2)/denom;
    % %     figure;plot(s, w', '-o'); title(strcat('lambda2 = ',num2str(opt.lambda2)));

end




% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 
% Visual samples in embedding space Y
% Required input:
%           .Y [nSmp 1]:    samples in embedding space (1-dimension)
%           .gnd [nSmp 1]:  class labels (start from 1)
%           .titlename:     figure title
% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 

function visualProj(Y, gnd, titlename)

    Y1 = Y(gnd == 1);
    Y2 = Y(gnd ~= 1);
    figure;     plot(Y1,ones(1,length(Y1)), 'bo');
    hold on;    plot(Y2,ones(1,length(Y2)), 'r+');
    title(titlename);

end




% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 
% Compute the center coordinate of each class in embedding space Y
% Required input m:
%           .Y [nSmp d]:    samples in embedding space (d-dimension)
%           .gnd [nSmp 1]:  class labels (start from 1)
% Output m:
%           .classCenter [nClass d]: coordinates of each class center
% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 

function m = classCenter(m)

    nClass = length(unique(m.gnd));
    m.classCenter = zeros(nClass, size(m.Y,2));
    for i = 1:nClass,
        m.classCenter(i,:) = mean(m.Y(m.gnd==i,:));
    end

end




% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 
% Classification with K-means (interface)
% Required input:
%       + data: sturct data type
%           .X [n]
%           .gnd
%       + iFold: 
% Output m:
%           .classCenter
% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 

function smodel = accuracyKmeans(data, iFold)

    smodel.predictedLabel = zeros(size(data.X,1),1);

    % split dataset to test and train data
    testIdx = (data.testSet == iFold);
    testData.X = data.X(testIdx,:);
    testData.gnd = data.gnd(testIdx,1);

    % classifiy test samples
    outlabel = binaryKmeans(data.classCenter, testData.X);
    % compute accuracy
    smodel.accuCV = sum(outlabel==testData.gnd)/sum(testIdx);
    % export predicted labels
    smodel.predictedLabel(find(data.testSet==iFold)) = outlabel;
    smodel.predictedLabel = [smodel.predictedLabel data.gnd data.testSet];

end




% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 
% Classification with K-means (binary case)
% Required input:
%       + classCenter [2 d]:    center coordinate of each class
%       + testX [nSmp d]:       test samples in embedding space
% Output:
%       + outlabel [nSmp 1]:    predicted labels (1 or 2)
% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - - 

function outlabel = binaryKmeans(classCenter, testX)
    
    [nSmp nDim] = size(testX);
    % distance from each sample to class centers
    classDist = sum((testX-repmat(classCenter(1,:),nSmp,1)).^2, 2) -...
                sum((testX-repmat(classCenter(2,:),nSmp,1)).^2, 2);
    % predict by the closest center
    outlabel = ones(nSmp, 1);
    outlabel(classDist > 0) = 2;

end




% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - -
% Evaluate SVM accuracy
% Required input:
%     + data: struct data type
%           .X [nFea nSmp]:     dataset X
%           .gnd [1 nSmp]:      class labels
%           .testSet [1 nSmp]:  sample indices for different CrossValidation
%     + idxFold: (optional)     if provided, compute accuracy for that fold only
% Output:
%     + model: struct data type
%           .accuCV:            X in embedding space
%           .Z:                 X in transformed space
%           .beta:              coeff to form each dim of transformed space
%           .predictedLabel:    predicted labels with additional cols of 
%                               groundtruth & CV-fold index
% - - - - - - - - - - - - - - - - - - - - - - - - - -  - - - - - - - - - - - - -

function m = accuracySVM(data, idxFold)

    if ~isfield(data,'testSet')
        error('This function requires cross validation indexing!');
    else
        bgnFold = 1;
        endFold = length(unique(data.testSet));
    end

    % If idxFold provided, perform a single fold train-test
    if (nargin==2), 
        bgnFold = idxFold;
        endFold = idxFold;
    end

    [nSmp, nFea] = size(data.X);
    if (bgnFold~=endFold)
        m.accuCV = zeros(1,endFold);
    end
    m.predictedLabel = zeros(1,nSmp);

    for iFold = bgnFold:endFold
    
        % split dataset to test and train data
        testIdx = (data.testSet==iFold);    trainIdx = ~testIdx;
        trainData.X = data.X(trainIdx,:);   trainData.gnd = data.gnd(trainIdx,1);
        testData.X = data.X(testIdx,:);     testData.gnd = data.gnd(testIdx,1);
        
        % train and test SVM
        outlabel = multisvm(trainData.X, trainData.gnd, testData.X);
        
        % compute accuracy
        if (bgnFold~=endFold)
            m.accuCV(iFold) = sum(outlabel==testData.gnd)/sum(testIdx);
        else
            m.accuCV = sum(outlabel==testData.gnd)/sum(testIdx);
        end

        % export predicted label 
        m.predictedLabel(find(data.testSet==iFold)) = outlabel;
    
    end
    
    m.predictedLabel = [m.predictedLabel' data.gnd data.testSet];

end
