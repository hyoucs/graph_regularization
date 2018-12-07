function linearSVM

f = load('data/maize_Data.mat');
dataAll = f.dataAll;
fout = fopen('record.txt', 'a+');

iData=2;
nFold = 5;

data = dataAll{iData};
%%nPosClass= length(find(data.gnd==1));
%%negClassIdx= find(data.gnd==0);

%%data.X(:,negClassIdx(nPosClass+1:end))=[];
%%data.gnd(:,negClassIdx(nPosClass+1:end))=[];
%%data.testSet(:,negClassIdx(nPosClass+1:end))=[];



%%rndIdx = randperm(length(negClassIdx), length(negClassIdx)-nPosClass);

%%data.X(:,negClassIdx(rndIdx)) = [];
%%data.gnd(:,negClassIdx(rndIdx)) = [];
%%data.testSet(:,negClassIdx(rndIdx)) =[];
CVO = cvpartition(data.gnd,'k',nFold); %stratified k-fold CV
for j=1:nFold
    teIdx = CVO.test(j);
    data.testSet (find(teIdx==1))=j;
    testArray{j}=teIdx';
end
nFold= length(unique(data.testSet));
figure; hist1D(data.gnd);
title(strcat(data.name,' percentage: ',num2str(sum(data.gnd)/length(data.gnd),'%4.2f')));
figure;
for j=1:nFold
    teIdx = (data.testSet == j);
    trIdx = ~teIdx;
    subplot(nFold,2,2*(j-1)+1);hist1D(data.gnd(teIdx));
    title(strcat(data.name,' percentage: ',...
        num2str(sum(data.gnd(teIdx))/length(data.gnd(teIdx)),'%4.2f')));
    subplot(nFold,2,2*j);hist1D(data.gnd(trIdx));
    title(strcat(data.name,' percentage: ',...
        num2str(sum(data.gnd(trIdx))/length(data.gnd(trIdx)),'%4.2f')));
end

dataDownSamples = data;
fprintf(fout, 'dataset: %s \n', data.name);
disp(data.name);

%%%%%%%%%%%% linear SVM on all features %%%%%%%%%%%%
data = dataDownSamples;
accuCV = lsvm(data);
disp(sprintf('all features: %5.2f| %4.2f', mean(accuCV),std(accuCV)));
fprintf(fout, 'all features: %5.2f| %4.2f \n', [mean(accuCV),std(accuCV)]);

 %%%%%%%%%%%% linear SVM on groundtruth features %%%%%%%%%%%%
 data=dataDownSamples;
 gndIdx = find(data.gndFea~=0);
 data.X = data.X(gndIdx, :);
 accuCV = lsvm(data);
 disp(sprintf('gnd features: %5.2f| %4.2f', mean(accuCV),std(accuCV)));
 fprintf(fout, 'rnd features: %5.2f| %4.2f \n', [mean(accuCV),std(accuCV)]);

%%%%%%%%%%%% linear SVM on randomized features %%%%%%%%%%%%
data=dataDownSamples;
gndIdx = find(data.gndFea~=0);
accuMX = [];
accuLT = [];
for j = 1:100,
    rndIdx = randperm(length(data.gndFea), length(gndIdx));
    rndata = data;
    rndata.X = rndata.X(rndIdx, :);
    accuCV = lsvm(rndata);
    accuMX = [accuMX, accuCV'];
    accuLT = [accuLT, mean(accuCV)];
    disp(sprintf('rnd features: %5.2f| %4.2f', mean(accuCV),std(accuCV)));
end
fprintf(fout, 'rnd features: %5.2f| %4.2f \n', [mean(accuLT)]);
disp(sprintf('average rnd features: %5.2f', mean(accuLT)));
hist(accuLT, 10);

fprintf(fout, '\n\n');


fclose(fout);

    function accuCV = lsvm(data)
        nFold = length(unique(data.testSet));
        accuCV = zeros(1,nFold);
        for iFold = 1:nFold %-- nFold CV
            testIdx = (data.testSet == iFold);
            trainIdx = ~testIdx;
            trainData.X = data.X(:,trainIdx);
            trainData.gnd = data.gnd(:,trainIdx);
            if isfield(data,'gndFea')
                trainData.gndFea = data.gndFea;
            end
            testData.X = data.X(:,testIdx);
            testData.gnd = data.gnd(:,testIdx);
            if (length(unique(data.gnd))==2)    %-- 2-class: use SVM
                svmStruct = svmtrain(trainData.X',trainData.gnd);
                outlabel = svmclassify(svmStruct,testData.X');
                accuCV(iFold) = sum(grp2idx(outlabel)==grp2idx(testData.gnd))/sum(testIdx);
            end
        end %-- nFold CV
    end

end
