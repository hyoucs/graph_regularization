function smodel = processOneFile(data,alphaArray,kNNArray)
%-- nFold cross validation over a single dataset
%-- cv_models ={m1 ... m10} for 10-fold CV
%-- param_model = {k alpha nfold_Models}
%-- mArray{end+1} = param_model
%-- smodel = mArray;


data.X = Xnorm(data.X);  
smodel={};

nAlpha = length(alphaArray);
nKnn = length(kNNArray);

nFold = length(unique(data.testSet));

options = [];
options.dDim=length(unique(data.gnd));
% options.dDim=5;

disp('kNN | alp |  AUC| accur| std');

for iKnn=1:nKnn  
    options.k = kNNArray(iKnn);
    %-- for each kNN/alpha, perform n-fold CV
    for iAlpha=1:nAlpha             
        options.alpha = alphaArray(iAlpha);
        v = [];

        %-- find U on entire DS
        fullModel = SNL(data,options);
        feaCoeff = max(abs(fullModel.U),[],2); 
        [dmp idx]= sort(feaCoeff,'descend');
        v.rankFea = [idx dmp]'; 
        
        %-- using default 50 topFea for transformed data
        nTopFea=50; 
        fullModel.U(idx(nTopFea+1:end),:)=0;
        data.Y = fullModel.U'*data.X;       %-- data for nFold-CV with SVM
        clear feaCoeff dmp idx

        v.kNN = kNNArray(iKnn);
        v.alpha = alphaArray(iAlpha);
        v.U = fullModel.U;
        v.eVal = fullModel.eVal;
        v.AUC = 0;
        if isfield(data,'gndFea') %-- compared to ground truth features
            if (isempty(v.eVal))
                error('Unable to test further on alpha and kNN');
                return
            end
            v.dROC = fullModel.dROC;
            v.AUC = prec_rec(fullModel.dROC(:,1),fullModel.dROC(:,2)); close;
        end
        
        %-- accuracy via cross validation
        v.accuCV = zeros(1,nFold);
        for iFold = 1:nFold %-- nFold CV
%             disp(sprintf('kNN %2d alpha %3.1f  %d-fold CV at: %d',...
%                 kNNArray(iKnn),alphaArray(iAlpha), nFold,iFold));
            %-- partition data
            testIdx = (data.testSet == iFold); 
            trainIdx = ~testIdx;
            
            trainData.X = data.X(:,trainIdx);
            trainData.Y = data.Y(:,trainIdx);
            trainData.W = data.W;
            trainData.gnd = data.gnd(trainIdx);
            if isfield(data,'gndFea')
                trainData.gndFea = data.gndFea;
            end
            
            testData.X = data.X(:,testIdx);
            testData.Y = data.Y(:,testIdx);
            testData.gnd = data.gnd(testIdx);


            if (length(unique(data.gnd))==2)    %-- 2-class: use SVM 
                svmStruct = svmtrain(trainData.Y',trainData.gnd);
                outlabel = svmclassify(svmStruct,testData.Y');
                v.accuCV(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);            
            else                                %-- >2-class: use k-means
                v.accuCV(iFold) = SNLpredict(testData,fullModel); 
            end
        end %-- nFold CV
        smodel = [smodel v];
        disp(sprintf('%3d |%5.1f| %4.2f| %5.2f| %4.2f',...
            v.kNN,v.alpha,v.AUC, mean(v.accuCV),std(v.accuCV)));
    end %-- Alpha loop
end %-- kNN loop 

end