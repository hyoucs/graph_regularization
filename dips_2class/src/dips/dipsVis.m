% visualize embedding space



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%				NON-LINEAR EMBEDDING
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% - - - - - - Load data
data_file = '../../../autism_private/data/correlation/cc200-fg-wcorr-whole-session-pc/dips_data.mat';
load(data_file);

meta_file = '../../../autism_private/data/nonlinear_output2.mat';
load(meta_file);

[pathstr,name,ext] = fileparts(data_file);
pw_file = [pathstr,'/dips_pwdist_',opt.type,ext];
load(pw_file);
data.A = pw.A;

nFold = length(unique(data.testSet));
idxFold = 2;

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

tr_1 = intersect(find(trainIdx==1), find(data.gnd==1));
tr_2 = intersect(find(trainIdx==1), find(data.gnd==2));
te_1 = intersect(find(testIdx==1), find(data.gnd==1));
te_2 = intersect(find(testIdx==1), find(data.gnd==2));
fine_gnd = zeros(size(data.gnd));
fine_gnd(tr_1) = 1; fine_gnd(tr_2) = 2;
fine_gnd(te_1) = 3; fine_gnd(te_2) = 4;

for idxl1 = 1:size(m2,1)
	for idxl2 = 1:size(m2,2)
		Xproj = data.X*m2{idxl1}{idxl2}.UAll;
		
		figure; 
		subplot(1,2,1); 
		gscatter(m.Y(:,1),m.Y(:,2),trainData.gnd);
		title('embeded');
		subplot(1,2,2);
		gscatter(Xproj(trainIdx,1),Xproj(trainIdx,2),trainData.gnd);
		title('approximated');

		figure; 
		subplot(1,2,1);
		gscatter(Xproj(trainIdx,1),Xproj(trainIdx,2),trainData.gnd);
		title('train approximated');
		subplot(1,2,2);
		gscatter(Xproj(:,1),Xproj(:,2),fine_gnd);
		title('test approximated');		
	end
end



% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 
%				LINEAR EMBEDDING
% - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - 

% - - - - - - Load data
data_file = '../../../autism_private/data/correlation/cc200-fg-wcorr-whole-session-pc/dips_data.mat';
load(data_file);	% data

eigs_file = '../../../autism_private/data/correlation/cc200-fg-wcorr-whole-session-pc/dips_linear_eig.mat';
load(eigs_file);	% x, d

meta_file = '../../../autism_private/data/meta/linear_opt.mat';
load(meta_file);	% opt

[pathstr,name,ext] = fileparts(data_file);
pw_file = [pathstr,'/dips_pwdist_',opt.type,ext];
load(pw_file);
data.A = pw.A;

nFold = length(unique(data.testSet));
idxFold = 1;

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

tr_1 = intersect(find(trainIdx==1), find(data.gnd==1));
tr_2 = intersect(find(trainIdx==1), find(data.gnd==2));
te_1 = intersect(find(testIdx==1), find(data.gnd==1));
te_2 = intersect(find(testIdx==1), find(data.gnd==2));
fine_gnd = zeros(size(data.gnd));
fine_gnd(tr_1) = 1; fine_gnd(tr_2) = 2;
fine_gnd(te_1) = 3; fine_gnd(te_2) = 4;


Xproj = data.X*d(:,1:2);

figure; 
subplot(1,2,1);
gscatter(Xproj(trainIdx,1),Xproj(trainIdx,2),trainData.gnd);
title('train approximated');
subplot(1,2,2);
gscatter(Xproj(:,1),Xproj(:,2),fine_gnd);
title('test approximated');		

