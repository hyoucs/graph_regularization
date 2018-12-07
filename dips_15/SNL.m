function model = SNL(data,options)
%------------------------------------------------------------------------%
% To find discriminative subnetworks to classify network samples 
%-- let u be eVec and a be eVal: 
%-- Min: u'[X(Lw-Lb)X' + alpha*C]u + beta*abs(u)
%-- Min: u'[X(Lw+Wb)X' + alpha*C]u + beta*abs(u)

%-- Max: u'[X(Lb-Lw)X' - alpha*C]U - beta*abs(u)
% Input:
%     + data: struct data type 
%           .X:   dataset X, size of [nFea nSmp]
%           .gnd: global network states, size of [1 nSmp]
%           .gndFea: global network states, size of [1 nSmp]
%           .W:   size of [nFea nFea], topology sharing among
%                 all networks xi's. wij is the frequency of edge ij in all networks
%     + options: structure 
%           .alpha: L2 tradeoff for smoothness topology, def.=0.1
%           .beta:  L1 tradeoff for sparness
%           .k:     number of nearest neighbors,    def.=5
%           .dDim:  number of leading eVecs kept for new data projection,
%                   def.= length(unique(gnd))
%           .nTopFea: top ranked features kept in eVecs, def.=All
% Output:
%     + model: struct data type 
%           .u: converged coefficients 
% Example:
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

tic
data = c2s(data);

Class = unique(data.gnd);
nClass = length(Class);

[nFea,nSmp] = size(data.X);

if (~exist('options','var'))
   options = [];
end

k = 5;
if isfield(options,'k'),     
    k = options.k; 
end

alpha = 0.1;
if isfield(options,'alpha'),    
    alpha = options.alpha; 
end

beta = 0.1;
if isfield(options,'beta'),    
    beta = options.beta; 
end

dDim=nClass;
if isfield(options,'dDim'),    
    dDim = options.dDim; 
end

nTopFea=nFea;
if isfield(options,'nTopFea'),    
    nTopFea = options.nTopFea; 
end

model.alpha = alpha;
model.k = k;
model.nTopFea = nTopFea;
model.kNN = cell(1,nSmp);     %-- cell for flexible # kNNs
model.U = zeros(nFea,dDim);   %-- store leading eVecs

%-- use cosine for kNN similarity
opt = [];
opt.bCosine = 1;
[dmp1 dmp2 A] = PWdistance(data.X,data.X,opt);
[dmp idx] = sort(A,2,'descend'); %-- 2: order in each row
clear dmp dmp1 dmp2

idx(:,[1,k+2:end]) = []; %-- note: 1st col-> self-node largest cosine
kNN = zeros(size(A));
for i = 1 : nSmp
    kNN(i,idx(i,:)) = 1;
end

A = A.*kNN;
A = max(A,A'); 

%-- build 2 meta-graphs: ML(within) and CL(btw) Affinity Matrices 

Aw = zeros(nSmp,nSmp);
Aw = sparse(Aw);
for i = 1 : nClass
    classIdx = find(data.gnd==Class(i));
    Aw(classIdx,classIdx) = 1;
end
Ab = ones(nSmp,nSmp)-Aw;

%-- Aw and Ab
Aw = Aw.*A; 
Ab = Ab.*A; 
Aw = max(Aw,Aw'); 
Ab = max(Ab,Ab'); 

%-- Lb= Db-Ab
Db = full(sum(Ab,2));
Lb = -Ab; clear Ab
for j = 1:size(Lb,1)
	Lb(j,j) = Db(j) + Lb(j,j) ;
end
clear Db

%-- Lw= Dw-Aw
Dw = full(sum(Aw,2));
Lw = -Aw; %clear Aw
for j = 1 : size(Lw,1)
	Lw(j,j) = Dw(j) + Lw(j,j) ;
end
Dw = diag(Dw);

%-- construct graph-constrainted Laplacian C=Lc=Dc-Wc
Wc = max(data.W,data.W'); 
Dc = full(sum(Wc,2));
C = -Wc; clear Wc
for j = 1 : size(C,1)
	C(j,j) = Dc(j) + C(j,j);
end

% %--- not much diff on syn data
% % C = diag(1./Dc)*C;%-- random walk Lalp
% C = diag(1./sqrt(Dc))*C*diag(1./sqrt(Dc));%-- symm Lalp
% C = max(C,C');
% %---
clear Dc


% % L = data.X*((beta/(1-beta))*Lb + Aw)*data.X'- alpha*C;   

%-- MAX case
%tildeL=data.X*(Lb - Lw)*data.X'; 
tildeL = data.X*(Lb + Aw)*data.X';
L = tildeL - alpha*C;  

% %-- MIN case (same results with eVal = diag(-eVal);)
% tildeL=data.X*(Lw-Lb)*data.X';
% L =  tildeL + alpha*C;  


L = max(L,L');
clear Lb Aw C 

if (nFea>=nSmp)    %-- singluarity
    Y= data.X*(Dw.^(.5)); %-- YY' = X(Dw^.5)(Dw^.5)X'
    [U,S,V]=SVD_CutOff(Y,SVDpct);
    invD=U*inv(S.^2)*U';
    invD = max(invD,invD');
    
%   issym = @(x) isequal(x,x.'); 
%   issym(invD*L) %-- not symm though both L and invD are!
%   length(find(invD*L<0))

    [eVec, eVal] = eigs(invD*L,dDim);   %-- solve D^-1Lu=a*u
else
    D = data.X*Dw*data.X';   
    D = max(D,D');
    [eVec, eVal] = eigs(L,D,dDim); %-- solve Lu=aDu
end

 
eVal = diag(eVal);
% eIdx = find(eVal < 1e-10);
% eVal(eIdx) = [];
% eVec(:,eIdx) = [];

dDim=min(dDim,length(eVal));
U=eVec(:,1:dDim); 
%-- Select nTopFea as final nodes 
feaCoeff = max(abs(U),[],2); 
[dmp idx]= sort(feaCoeff,'descend'); 
U(idx(nTopFea+1:end),:)=0;

model.U=U; 
model.eVal=eVal;
% model.feaCoeff=feaCoeff; 
% model.feaIdx=idx; 
model.dDim=dDim;
model.elptime=toc;


%-- compute dROC if gndFea is provided
if isfield(data,'gndFea') 
    if  (length(idx)>0),
        model.dROC=[feaCoeff(idx) data.gndFea(idx)'];   
    else
        model.dROC=[];        
    end
end

%-- Class means for cross-validation testing
Y=model.U'*data.X; %-- in transformed space
% figure;PlotX(model.Y,gnd,'','Samples in transformed space','','');

ClassCenter=zeros(dDim,nClass);
for i=1:nClass
    iClass = Y(:,data.gnd==Class(i));
    ClassCenter(:,i) = mean(iClass,2);
end
model.ClassCenter = ClassCenter;
model.ClassLabel = Class;
model.gnd = data.gnd;
model.algo = 'SNL';

end