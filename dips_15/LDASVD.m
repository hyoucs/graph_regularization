function model = LDASVD(data,opt)
%------------------------------------------------------------------------%
% Perform LDA for cases D>>N using pseudo inverse
% Input:
%     + data: 
%         .X [nFea nSmp]: data set
%         .gnd [1 nSmp]: class labels
%     + opt:
%         .nDim: new subspace's dim, default == nClass
%         .pct: % of singular values kept (btw (0,1])
%             =1: Sw is assumed full rank, nonsingluar, but SVD can still
%             be used, this case U=V
%             <1: select top singluars that occupy pct total singluars.
%             SVD is used to approximate Sw. 
%         .verbose: 1/0 for display results, default 0-not display
% Output:
%      + model: includes Sb,Sw,W,... for direction and projection
% Note:
%    + Objective: Sw can be not full-rank by checking Sw=U*S*V'. Thus, obj
%    is to approx inv(Sw)= V'*inv(S)*U. eig(inv_Sw*Sb) does not affected!     
%    + Codes adapted from stprtool

% (c) 2012 Quang-Anh Dang - Aarhus University 
%     2014-09-04: udpated at UCSB
%------------------------------------------------------------------------%


[nFea nSmp] = size(data.X);
sClass = unique(data.gnd);
nClass = length(sClass);


nDim = nClass;
if isfield(opt,'nDim'),
    nDim = opt.nDim;
end

pct = 1;
if isfield(opt,'pct'),
    pct = opt.pct;
end


verbose = 0;
if isfield(opt,'verbose'),
    verbose = opt.verbose;
end



Sw = zeros(nFea,nFea);
Sb = Sw;

meanX = mean(data.X,2);

for i = 1:nClass,
  idx_i = find(data.gnd==sClass(i));
  ni = length(idx_i);
  X_i = data.X(:,idx_i);
  meanXi = mean(X_i,2);
  
  if size(X_i,2)>1 %-- compute cov.mat. for class with at least > 1 subject
    Sw = Sw + cov(X_i',1);
  end
    
  Sb = Sb + ni*(meanXi-meanX)*(meanXi-meanX)';%-- DQA: weighted by class size!!!!
%   Sb = Sb + (mean_Xi-mean_X)*(mean_Xi-mean_X)';%-- DQA: no weight
end



if (pct == 1)    % Sw is full rank (nonsingluar)
    %-- NOTE: no difference if SVD is used for FULL_RANK. This case U=V   
    [U,S,V]=svd(Sw);
    model.S = diag(S);
    
    if ~isempty(find(diag(S)<0.0000001))
        error('Within matrix is singular! Better to set pct<1');
    end
    
    inv_Sw=V*inv(S)*U'; %-- pseudo inverse 
    [V,D]=eig(inv_Sw*Sb);
    S = diag(S);
   
    %-- or simply perform eig func:
    %[V,D]=eig(inv(Sw)*Sb); %-- usual case    

else % approx Sw: keeping pct of total singulars    
    [U,S,V]=svd(Sw); %-- then Sw=U*S*V' and inv(Sw)=V*inv(S)*U'
    S = diag(S);
    % select top singluars
    for i=1:length(S)
        sum(S(1:i))/sum(S);
        if(sum(S(1:i))/sum(S)>=pct) break; end
    end   
    
    model.S = S(1:i);
    S = diag(S(1:i));
    inv_Sw=V(:,1:i)*inv(S)*U(:,1:i)'; %-- pseudo inverse 
    [V,D]=eig(inv_Sw*Sb);
end

[D,inx] = sort(diag(D),1,'descend');

% take nDim biggest eigenvectors
model.V = V(:,inx(1:nDim)); %-- set of leading eVecs

% translation
model.b = -model.V'*meanX; %-- projection of mean on W
model.X = model.V'*data.X + model.b(:)*ones(1,nSmp);
% model.X = -model.V'*data.X; %-- similar effect

model.Sw = Sw;
model.Sb = Sb;
model.meanX = meanX;

model.fun = 'LDA with pseudo inverse';


if verbose
    figure;
    subplot(1,3,1); PlotX(data.X,data.gnd,'','Original data with groundtruth','',''); grid on;
    subplot(1,3,2); plot(model.S,'b+'), title ('Eigen Values on SVD'); grid on; hold on; xlim([0 length(S)+1]);
    subplot(1,3,3); PlotX(model.X,data.gnd,'','data in spectral space, from 2nd eVec','',''); grid on;
end


return;
% EOF