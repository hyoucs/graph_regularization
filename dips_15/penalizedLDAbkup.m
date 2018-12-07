function model = penalizedLDAbkup(data,options,varargin)
%------------------------------------------------------------------------%
% Implement penalizedLDA following Witten and Irina paper
%-- Max: v'Bv - beta*abs(v) s.t. v'Wv <=1 
%-- 2 steps: 
%--     u = B^{.5}*v 
%--     v= q/(q'*W*q) where q is found from 2*B^{.5}*u - 2W*q - softfunc
%
% Input:
%     + data: struct data type
%           .X:   dataset X, size of [nFea nSmp] .gnd: global network
%           states, size of [1 nSmp] .gndFea: global network states, size
%           of [1 nSmp] .W:   size of [nFea nFea], topology sharing among
%                 all networks xi's. wij is the frequency of edge ij in all
%                 networks
%     + options: structure
%           .lambda1: L1 tradeoff for sparness .lambda2: L2 tradeoff for
%           smoothness topology, def.=0.1 .k:     number of nearest
%           neighbors,    def.=5 .dDim:  number of leading eVecs kept for
%           new data projection,
%                   def.= length(unique(gnd))
%           .nTopFea: top ranked features kept in eVecs, def.=All
% Output:
%     + model: struct data type
%           .u: converged coefficients
% Example: (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

tic
data = c2s(data);

Class  = unique(data.gnd);
nClass = length(Class);

[nFea,nSmp] = size(data.X);

if (~exist('options','var'))
    options = [];
end

lambda1 = 0.1;
if isfield(options,'lambda1'),
    lambda1 = options.lambda1;
end

dDim=nClass-1;
if isfield(options,'dDim'),
    dDim = options.dDim;
end



%% -- Compute B and W by LDA
data.X = data.X'; %-- transpose X for following computation
% Make sure data is zero mean
mapping.mean = mean(data.X, 1);
data.X = bsxfun(@minus, data.X, mapping.mean);

% Make sure labels are nice
[Class, bar, labels] = unique(data.gnd);
nClass = length(Class);

% Intialize Sw
Sw = zeros(nFea, nFea);

% Compute total covariance matrix
St = cov(data.X);

% Sum over classes
for i=1:nClass
    cur_X = data.X(labels == i,:);
    C = cov(cur_X);
    p = size(cur_X, 1) / (length(labels) - 1);
    Sw = Sw + (p * C);
end

% Compute between class scatter
Sb       = St - Sw;
Sb(isnan(Sb)) = 0; Sw(isnan(Sw)) = 0;
Sb(isinf(Sb)) = 0; Sw(isinf(Sw)) = 0;

W = Sw; B = Sb; clear Sb Sw St cur_X
data.X = data.X';

%% -- end computing B, W and initial theta by LDA


[M, lambda] = eig(B, W);
% Sort eigenvalues and eigenvectors in descending order
lambda(isnan(lambda)) = 0;
[lambda, ind] = sort(diag(lambda), 'descend');
M = M(:,ind(1:min([dDim size(M, 2)])));



%%
theta = M(:,1); %-- initialize v as 1st eVec found from LDA
% theta = rand(nFea,1); %-- unstable results!

[maxIter,verbose,optTol,zeroThreshold] = process_options([],'maxIter',1000,'verbose',2,'optTol',1e-5,'zeroThreshold',1e-4);

verbose=1;
if verbose==2       %Start the log
    fprintf('%10s %10s %15s %15s %15s\n','iter','shoots','abs(bi)','abs(bi-bj)','f(w)');
    theta_log(:,niter+1) = theta;
end

bFullW=0;    %-- full W, not run perfectly!
if (~bFullW) %-- diag W, thus a vector
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
