function [m] = accuracySVM(data,idxFold)
%------------------------------------------------------------------------%
% Evaluate SVM accuracy
% Input:
%     + data: struct data type
%           .X [nFea nSmp]: dataset X
%           .gnd [1 nSmp]: class labels
%           .testSet [1 nSmp]: sample indices for different CrossValidation
%     + idxFold: optional
%           . if provided, compute accuracy for that fold only
% Output:
%     + model: struct data type
%           .accuCV: X in embedding space
%           .Z: X in transformed space
%           .beta: coeff to form each dim of transformed space
%           .predictedLabel: predicted labels with additional cols of groundtruth & CV-fold index
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%


if ~isfield(data,'testSet')
    error('This function requires cross validation indexing!');
else
    bgnFold = 1;
    endFold = length(unique(data.testSet));
end


if (nargin==2), 
    bgnFold = idxFold;
    endFold = idxFold;
end

[nFea nSmp] = size(data.X);

m.predictedLabel = zeros(1,nSmp);

if (bgnFold~=endFold)
    m.accuCV = zeros(1,endFold);
end

for iFold = bgnFold:endFold
    % split dataset to test and train data
    testIdx = (data.testSet == iFold);
    trainIdx = ~testIdx;
    
    trainData.X = data.X(:,trainIdx);
    trainData.gnd = data.gnd(trainIdx);
    
    testData.X = data.X(:,testIdx);
    testData.gnd = data.gnd(testIdx);
    
    % train and test SVM
    outlabel = multisvm(trainData.X',trainData.gnd',testData.X');
    if (bgnFold~=endFold)
        m.accuCV(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
    else
        m.accuCV = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
    end
    m.predictedLabel(find(data.testSet == iFold)) = outlabel;
end
m.predictedLabel = [m.predictedLabel' data.gnd' data.testSet'];

end
