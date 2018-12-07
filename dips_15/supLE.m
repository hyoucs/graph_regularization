function model = supLE(data,options,varargin)
%------------------------------------------------------------------------%
%-- Supervised Laplacian Embedding techniques, both linear and non-linear projections 
%-- Non-linear: Max y'(Lb + Aw)y s.t. y'Dwy=1  
%-- Linear:     a'X(Lb + Aw)X'a s.t. a'XDwX'a=1  (might add a'Ca as network topo constraint!)
% 
%-- Note: Different sizes of eigVecs y (1*nSmp) and a (1*nFea)
%     Non-linear: entries in y are mapping points
%     Linear: entries in a are linear coeficients 

% To find discriminative subnetworks to classify network samples
% 2 steps:
% seek y s.t. MAX y'(Lb+Aw)y s.t. y'Dwy=1 or equivalently
%               z'(Dw^{-.5}(Lb+Aw)Dw^{-.5})z = z'Lz s.t. z'z=1
%  view L=XX' and decompose X=USV', similar to sparsePCA
%  find sparse a s.t. min ||Xa -y ||^2  + \alpha a'Ca + \beta |a|
%  where y' = \sigma*v' by svd


% Input:
%     + data: struct data type
%           .X [nFea nSmp]:   dataset X
%           .gnd [1 nSmp]: global network states
%           .gndFea [1 nSmp]: global network states
%           .W [nFea nFea]:topology sharing among samples
%                 Wij is the frequency of edge ij in all networks
%     + options: structure
%           .bLinear: =1 for linear case
%           .k:     number of nearest neighbors,    def.=5
%           .alpha: tradeoff (alpha*Lb + (1-alpha)*Aw), def.=.5
%           .lambda2: L2 tradeoff for smoothness topology, def.=0.1
%           .dDim:  projected subspace dimension def.= nClass-1
% Output:
%     + model: struct data type
%           .eVecs: either y or a above. 
%           .eVals: either y or a above. 
%           .Y: samples presented in projected subspace (=y in non-linear,
%           = a'X in linear)
%           .Lb = Db-Ab
%           .Aw 
%           .Dw 
% Example:
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

tic
data = c2s(data);

Class  = unique(data.gnd);
nClass = length(Class);

[nFea,nSmp] = size(data.X);

if (~exist('options','var'))
    options = [];
end

k = 10;
if isfield(options,'k'),
    k = options.k;
end

bLinear = 1;
if isfield(options,'bLinear'),
    bLinear = options.bLinear;
end


alpha = 0.5;
if isfield(options,'alpha'),
    alpha = options.alpha;
end

lambda2 = 0.1;
if isfield(options,'lambda2'),
    lambda2 = options.lambda2;
end

dDim=nClass-1;
if isfield(options,'dDim'),
    dDim = options.dDim;
end

model.alpha = alpha;
model.k = k;

%-- use Euclidean for kNN similarity (or Pearson for expression profiles)
%-- Q: for gene data: is 0-corr/cos gene-pair more similar compared to neg-corr gene-pair?
opt = [];
opt.type='Eucl';
% opt.type='Cosine';
[A] = PWdistance(data.X,data.X,opt);
[dmp idx] = sort(A,2,'descend'); 

idx(:,[1,k+2:end]) = []; %-- remove 1st (selfnode sim) + last smallest sim
kNN = zeros(size(A));
for i = 1 : nSmp
    kNN(i,idx(i,:)) = 1;
end

A = A.*kNN;
A = max(A,A'); %-- (attention to corr/cosine similarity)

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

% figure;imshow(full(A),'InitialMagnification',2000);


%-- Compute Lb= Db-Ab
Db = full(sum(Ab,2));
Lb = -Ab; %clear Ab
for k = 1:size(Lb,1)
    Lb(k,k) = Db(k) + Lb(k,k) ;
end
clear Db

%-- Compute Lw= Dw-Aw
Dw = full(sum(Aw,2));
if (length(find(Dw==0))>0),
    error('Increase kNN to ensure no disconnected "within" nodes!');
end

Lw = -Aw; %clear Aw
for k = 1 : size(Lw,1)
    Lw(k,k) = Dw(k) + Lw(k,k);
end
Dw = diag(Dw);


% %-- construct graph-constrainted Laplacian C=Lc=Dc-Wc
% Wc = max(data.W,data.W');
% Dc = full(sum(Wc,2));
% C = -Wc; clear Wc
% for k = 1 : size(C,1)
%     C(k,k) = Dc(k) + C(k,k);
% end
%
% clear Dc

%-- MAX y'(Lb+Aw)y s.t. y'Dwy=1 or z'Lz s.t. z'z=1 with z=Dw^{-.5}y
%-- alpha: affects L1 selection,
%       better for small alpha ~ more emphasis on Aw
%       if balanced alpha, kNN should be large, say 40
% alpha=0.1;
% L=(alpha*Lb + (1-alpha)*Aw); %-- L=(Lb - Lw);  %-- returns same result
% L = max(L,L');
% diag(L)

% alpha=0.5; %-- Aw is much more important than Lb
% if bLinear %-- linear case: max a'X(Lb+Aw)X'a s.t. a'XDX'a=1
%    L = data.X*(alpha*Lb + (1-alpha)*Aw)*data.X';
%    Dw = data.X*Dw*data.X';
% else       %-- non-linear case: max a'(Lb+Aw)a s.t. a'Da=1
%     L = alpha*Lb + (1-alpha)*Aw;
% end



L = alpha*Lb + (1-alpha)*Aw;
% L = Lb + alpha*Aw;

% linear case (slower): max a'X(Lb+Aw)X'a s.t. a'XDX'a=1
if bLinear 
   L = data.X*L*data.X';
   Dw = data.X*Dw*data.X';
end



% %-- pseudo inv is worse than using eig
% [U,S,V]=SVD_CutOff(Dw,.99);
% invD = V*inv(S)*U';
% invD = max(invD,invD');
% [eVecs, eVals] = eigs(invD*L,dDim); 
% diag(eVals)
% Dw*invD

% %-- better for singular Dw, check algo behind eig and eigs
% [eVecs, eVals] = eig(L, Dw); % diag(eVals) % work well for non-linear but
% show for linear case on Brain data



[eVecs, eVals] = eig(L, Dw);
eVals(isnan(eVals)) = 0;
[eVals, ind] = sort(diag(eVals), 'descend');
eVecs = eVecs(:,ind(1:min([dDim size(eVecs, 2)])));
  


model.eVecs = eVecs;
model.eVals = eVals;
if bLinear 
    model.Y = eVecs'*data.X;
else       
    model.Y = eVecs';
end
model.Lb = Lb;
model.Aw = Aw;

model.B = L;
model.W = Dw;

% figure;
%     subplot(1,3,1); PlotX(data.X,data.gnd,'','Original space',''); grid on;
%     subplot(1,3,2); PlotX(model.Y,data.gnd,'','Embedded space','');
%     subplot(1,3,3); stem(abs(model.eVecs(:,1)));

end