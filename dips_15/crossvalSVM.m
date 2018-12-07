function [accSum accFull] = crossvalSVM(data,isGTfea)
%------------------------------------------------------------------------%
%-- perform SVM on n-folds and report accuracy
%-- data must follow std format
% if isGTfea: gndFea is provided, perform 3 tests, 1 for entire
% Input:
%     + data: struct data type (std format as that for SNL)
%           .X:   dataset X, size of [nFea nSmp]
%           .gnd: global network states, size of [1 nSmp]
%           .testSet: for cross validation
%     + isGTfea (optional): further test on GTFeature and Random features
%     (same size)
% Output:
%     + accSum: mean + stddev for each case
%     + accFull: accuracy of each individual run
% Example:
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%



%-- full dataset
acc = cvSVM(data);
accSum = [mean(acc) std(acc)];
accFull = acc;


if (exist('isGTfea','var'))

   %-- run on GT-feature
   dataGT = data;
   gtFeaIdx = find(data.gndFea==1);
   nGTfea = length(gtFeaIdx);
   
   dataGT.X = data.X(gtFeaIdx,:); 
   acc = cvSVM(dataGT);
   clear dataGT;
   
   accSum  = [accSum mean(acc) std(acc)];
   accFull = [accFull acc];
  
   
   %-- run on features randomly selected from non-GTfea 
   %-- (same size as GTfea)
   
   dataRD = data;
   rdFeaIdx = find(data.gndFea==0);
   acc=[];
   for iRand = 1:10 %-- number of random running
       posSwap = randperm(length(rdFeaIdx));
       rdFea = rdFeaIdx(posSwap(1:nGTfea));
     
       dataRD.X = data.X(rdFea,:);
       acc =  [acc cvSVM(dataRD)];
   end
   
   clear dataRD;
   
   accSum  = [accSum mean(acc) std(acc)];
   accFull = [accFull acc];

end


function [accu] = cvSVM(data)
%-- function run n-fold cv for each dataset

nFold = length(unique(data.testSet));
accu = zeros(1,nFold);

for iFold = 1:nFold %-- nFold CV
    %-- partition data
    testIdx = (data.testSet == iFold);
    trainIdx = ~testIdx;
    
    trainData.X = data.X(:,trainIdx);
    trainData.gnd = data.gnd(trainIdx);
    testData.X = data.X(:,testIdx);
    testData.gnd = data.gnd(testIdx);
    
    svmStruct = svmtrain(trainData.X',trainData.gnd);
    outlabel = svmclassify(svmStruct,testData.X');
    accu(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
   
end %-- nFold CV
