function model = RegSubnetNoConstraint(data,options,varargin)
%------------------------------------------------------------------------%
% To find discriminative subnetworks to classify network samples
% Not done as u = 0 leads to minimum due to lack constraints imposed on u!
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

Class  = unique(data.gnd);
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
Lw = -Aw; %clear Aw
for k = 1 : size(Lw,1)
    Lw(k,k) = Dw(k) + Lw(k,k);
end
Dw = diag(Dw);

%-- construct graph-constrainted Laplacian C=Lc=Dc-Wc
Wc = max(data.W,data.W');
Dc = full(sum(Wc,2));
C = -Wc; clear Wc
for k = 1 : size(C,1)
    C(k,k) = Dc(k) + C(k,k);
end

% %--- not much diff on syn data
% % C = diag(1./Dc)*C;%-- random walk Lalp
% C = diag(1./sqrt(Dc))*C*diag(1./sqrt(Dc));%-- symm Lalp
% C = max(C,C');
% %---
clear Dc


% % L = data.X*((beta/(1-beta))*Lb + Aw)*data.X'- alpha*C;

%-- from here!


%-- MIN case

% tildeL=data.X*(Lb - Aw)*data.X';
tildeL=data.X*(Lw + Ab)*data.X';
L = tildeL + 0.01*eye(nFea); %-- 0.01*C
L = max(L,L'); %-- check diagonals >0?
% diag(L)
clear Lb Aw C tildeL


[maxIter,verbose,optTol,zeroThreshold] = process_options(varargin,'maxIter',10000,'verbose',2,'optTol',1e-5,'zeroThreshold',1e-4);

niter = 0;
theta = rand(nFea,1);
% theta = (data.X*data.X' + beta*eye(nFea))\(data.X*data.gnd');
if verbose==2       %Start the log
    fprintf('%10s %10s %15s %15s %15s\n','iter','shoots','abs(bi)','abs(bi-bj)','f(w)');
    theta_log(:,niter+1) = theta;
end


while niter < maxIter
    theta_old = theta;
    
    for k = 1:nFea %-- update w.r.t. each feature
        
        ak = 2*L(k,k);
        ck = 2*(theta(k)*L(k,k) - theta'*L(:,k));
        
        if ck > beta
            theta(k) = (ck-beta)/ak;
        elseif ck < -beta
            theta(k) = (ck+beta)/ak;
        elseif abs(ck) <= beta
            theta(k) = 0;
        end        
        
    end
    
    niter = niter + 1;
    
    % Update the log
    if verbose==2
        fprintf('%10d %10d %15.4e %15.4e %15.4e\n',niter,niter*nFea,sum(abs(theta)),...
            sum(abs(theta - theta_old)),...
            sum((data.X'*theta - data.gnd').^2)+ beta*sum(abs(theta)));
        theta_log(:,niter+1) = theta;
    end
    % Check termination
    if sum(abs(theta-theta_old)) < optTol
        break;
    end
    
    
end
% if verbose
%     fprintf('#iterations: %d  -- Total Shoots: %d\n',niter-1,(niter-1)*nFea);
% end
model.u = theta;
model.theta_log = theta_log;
%-- end from here


end