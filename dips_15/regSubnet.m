function model = regSubnet(data,options,varargin)

%------------------------------------------------------------------------%
% To find discriminative subnetworks to classify network samples
%-- Max: v'[X(Lb-Lw)X' - alpha*C]v - beta*abs(v) s.t. v'*X*Dw*X'*v <=1
%-- Max: v'Bv - beta*abs(v) s.t. v'Wv <=1 
%-- 2 steps: u = B^{.5}*v
%--  v= q/(q'*Dw*q) where q is found from 2*B^{.5}*u - 2Dw q - softfunc
% 
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


lambda2 = 0.1;
if isfield(options,'lambda2'),
    lambda2 = options.lambda2;
end

lambda1 = 0.1;
if isfield(options,'lambda1'),
    lambda1 = options.lambda1;
end


bLinear = 1;
if isfield(options,'bLinear'),
    bLinear = options.bLinear;
end


lambda1 = 0.1;
if isfield(options,'lambda1'),
    lambda1 = options.lambda1;
end

lambda2 = 0.1;
if isfield(options,'lambda2'),
    lambda2 = options.lambda2;
end

dDim=nClass-1;
if isfield(options,'dDim'),
    dDim = options.dDim;
end

nTopFea=nFea;
if isfield(options,'nTopFea'),
    nTopFea = options.nTopFea;
end

model.lambda2 = lambda2;
model.k = k;
model.nTopFea = nTopFea;
model.kNN = cell(1,nSmp);     %-- cell for flexible # kNNs

%-- use Euclidean for kNN similarity (or Pearson for expression profiles)
%-- Q: for gene data: is 0-corr/cos gene-pair more similar compared to neg-corr gene-pair?
opt = [];
opt.bEucl = 1;
[dmp1 A dmp2 dmp3] = PWdistance(data.X,data.X,opt);

% opt.bGauss = 1; opt.bEucl=1;
% [A dmp1 dmp2 dmp3] = PWdistance(data.X,data.X,opt);

[dmp idx] = sort(A,2,'descend'); %-- 2: order in each row: 0-sim >> neg.sim??
clear dmp dmp1 dmp2 dmp3

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
% figure;imshow(full(Aw),'InitialMagnification',2000);
% figure;imshow(full(Ab),'InitialMagnification',2000);


%-- Lb= Db-Ab
Db = full(sum(Ab,2));
Lb = -Ab; %clear Ab
for k = 1:size(Lb,1)
    Lb(k,k) = Db(k) + Lb(k,k) ;
end
clear Db

%-- Lw= Dw-Aw
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

alpha=0.5; %-- Aw is much more important than Lb
if bLinear %-- linear case: max a'X(Lb+Aw)X'a s.t. a'XDX'a=1
   L = data.X*(alpha*Lb + (1-alpha)*Aw)*data.X';
   Dw = data.X*Dw*data.X';
else       %-- non-linear case: max a'(Lb+Aw)a s.t. a'Da=1
    L = alpha*Lb + (1-alpha)*Aw;
end


[eVec, eVal] = eigs(L,Dw,dDim);

model.eVec = eVec;
model.eVal  = eVal;


if bLinear 
    model.y = eVec'*data.X;
else       
    model.y = eVec';
end

B = L; W = Dw; 

%%
theta = eVec(:,1); %-- initialize v as 1st eVec found
% theta = rand(nFea,1); %-- unstable results!

[maxIter,verbose,optTol,zeroThreshold] = process_options([],'maxIter',1000,'verbose',2,'optTol',1e-5,'zeroThreshold',1e-4);

verbose=1;
if verbose==2       %Start the log
    fprintf('%10s %10s %15s %15s %15s\n','iter','shoots','abs(bi)','abs(bi-bj)','f(w)');
    theta_log(:,niter+1) = theta;
end

bFullW=0;
if (~bFullW) %-- simplify W as diag.mat and thus a vector
    W=diag(W); 
end

[U,S,V]=SVD_CutOff(B,1);
Bhalf=U*(S.^.5)*V';

niter = 0;
qold = theta;
qnew = zeros(nFea,1);
while niter < maxIter
    theta_old = theta;
    %-- 1st step compute u
    u = Bhalf*theta_old;
%     [lambda1/2 max(u)]
    
    %-- compute theta (or v) through q
    
    for k = 1:nFea %-- coordinate ascent q
        utB = u'*Bhalf(:,k); 
        if (bFullW) %-- full W case
            wk = W(k,k);
            utB = utB - qold'*W(:,k) + qold(k)*wk;
        else
            wk = W(k);
        end
        
        qnew(k)=softThreshold(utB,lambda1/2)/wk;
        
%         if utB > lambda1/2
%             q(k) = (utB - lambda1/2)/W(k);
%         elseif utB < -lambda1/2
%             q(k) = (utB + lambda1/2)/W(k);
%         elseif abs(utB) <= lambda1/2
%             q(k) = 0;
%         end
    end
    qold=qnew;
    
    %-- update theta based on q
    if ~any(qnew) %-- ALL entries are 0
        theta = zeros(nFea,1);
    else
        if (bFullW) %-- full W case
            theta = qnew/(sqrt(qnew'*W*qnew)); 
        else        
            theta = qnew/(sqrt(qnew'*(W.*qnew))); 
        end
    end
    
    niter = niter + 1;
    % Update the log
    if verbose==2
        fprintf('%10d %10d %15.4e %15.4e %15.4e\n',niter,niter*nFea,sum(abs(theta)),...
            sum(abs(theta - theta_old)),...
            sum((data.X'*theta - data.gnd').^2)+ lambda1*sum(abs(theta)));
        theta_log(:,niter+1) = theta;
    end
    % Check termination
    if sum(abs(theta-theta_old)) < optTol
        %         fprintf('diff: %f',sum(abs(theta-theta_old)));
        break;
    end
end

% if verbose
%     fprintf('#iterations: %d  -- Total Shoots:
%     %d\n',niter-1,(niter-1)*nFea);
% end

model.v = theta;
% model.theta_log = theta_log;


%-------------------------------------------------------------------------%

function [scal]=softThreshold(mat,lam)
    %-- compute soft func, both mat and lam are scalars
    scal = sign(mat)*max(abs(mat)-lam, 0);
end

end