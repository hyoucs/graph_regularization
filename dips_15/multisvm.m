function [result] = multisvm(TrainingSet,GroupTrain,TestSet)
%Models a given training set with a corresponding group vector and 
%classifies a given test set using an SVM classifier according to a 
%one vs. all relation. 


u=unique(GroupTrain);
numClasses=length(u);
result = zeros(length(TestSet(:,1)),1);

% replace svmtrain with fitcsvm, starting from Matlab 2018a
% replace 
%             VMStruct = svmtrain(X,Y);
%             label = svmclassify(SVMStruct,Xnew);
% with
%             SVMModel = fitcsvm(X,Y);
%             label = predict(SVMModel,Xnew);

models = {};

%build models
for k=1:numClasses
    %Vectorized statement that binarizes Group
    %where 1 is the current class and 0 is all other classes
    G1vAll=(GroupTrain==u(k));
    models{k} = fitcsvm(TrainingSet,G1vAll);
end

%classify test cases
for j=1:size(TestSet,1)
    for k=1:numClasses
        if(predict(models{k},TestSet(j,:))) 
            break;
        end
    end
    result(j) = k;
end