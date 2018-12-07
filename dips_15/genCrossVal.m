function [data]= genCrossVal(data,nFold)
%------------------------------------------------------------------------%
% Call this func to generate "nFold" stratified cross validation
% Input:
%     + data: struct data type
%           .X [nFea nSmp]: dataset X (no need?)
%           .gnd [1 nSmp]: class labels
%     + nFold: number of folds
% Output:
%     + data: with testSet feature for CV
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%


testArray=cell(1,nFold);
data.testSet=zeros(1,length(data.gnd));
figure; 
hist1D(data.gnd);
title(' class percentage overall dataset');
CVO = cvpartition(data.gnd,'k',nFold); %stratified k-fold CV
figure; 
for j=1:nFold
    teIdx = CVO.test(j);
    trIdx = ~teIdx;
    subplot(nFold,2,(j-1)*2+1);hist1D(data.gnd(teIdx));
    title('Test data');
    
    subplot(nFold,2,(j-1)*2+2);hist1D(data.gnd(trIdx));
    title('Train data');
    data.testSet (find(teIdx==1))=j;
    testArray{j}=teIdx';
end
end