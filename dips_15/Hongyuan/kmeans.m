function [m] = kmeans(data,idxFold)
%KMEANS Summary of this function goes here
%   Detailed explanation goes here

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
    
% %     % libsvm
%     model = svmtrain(trainData.gnd', trainData.X', '-t 0 -q -nu 0.1');
%     [outlabel, ~, ~] = svmpredict(testData.gnd', testData.X', model, '-q');
%     v.accuCV(iFold) = accuracy(1)*0.01;

%     % matlab svm
%     svmStruct = svmtrain(trainData.X',trainData.gnd);
%     outlabel = svmclassify(svmStruct,testData.X');
% 	  v.accuCV(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
    
    % train and test SVM
    outlabel = multisvm(trainData.X',trainData.gnd',testData.X');
    if (bgnFold~=endFold)
        m.accuCV(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
    else
        m.accuCV = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
    end
    m.predictedLabel(find(data.testSet == iFold)) = outlabel;
end


end

