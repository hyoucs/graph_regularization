function model = coordinateAscent(data,options,varargin)
%------------------------------------------------------------------------%
% Max [v'Bv - beta*abs(v)] s.t. v'Wv <=1 
%               with B, W (e.g., by LDA) and given initial beta
%-- 2 steps: 
%--     u = B^{.5}*v 
%--     v = q/(q'*W*q) where q is found from 2*B^{.5}*u - 2W*q - softfunc
%
% Input:
%     + data: struct data type
%           .X[nFea nSmp]: dataset X (no need?)
%           .gnd[1 nSmp]: global network states (no need?)
%           .B[nFea nFea]: between-class matrix 
%           .W[nFea nFea]: within-class covariance   
%     + options: structure
%           .lambda1:   L1 tradeoff for sparsness 
%           .lambda2:   L2 tradeoff for smoothness
%           .theta [nFea 1]: initial theta
%           .bFullW: boolean, using full or diag. W
%           .verbose: boolean, used to print out results
% Output:
%     + model: struct data type
%           .v: sparse vector of coeffiences
% Example: (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%


%%

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


theta = rand(nFea,1); %-- unstable results!
if isfield(options,'theta'),
    theta = options.theta; 
end

bFullW=0;    %-- full W, not run perfectly!
if isfield(options,'bFullW'),
    bFullW = options.bFullW; 
end


verbose = 0;
if isfield(options,'verbose'),
    verbose = options.verbose; 
end


dDim=nClass-1;
if isfield(options,'dDim'),
    dDim = options.dDim;
end

B = data.B;
W = data.W;


%[maxIter,verbose,optTol,zeroThreshold] = process_options([],'maxIter',1000,'verbose',2,'optTol',1e-5,'zeroThreshold',1e-4);

maxIter = 1000;
optTol = 1e-5;
zeroThreshold = 1e-4;

if verbose       %Start the log
    fprintf('%10s %10s %15s %15s %15s\n','iter','shoots','abs(bi)','abs(bi-bj)','f(w)');
    theta_log(:,niter+1) = theta;
end

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
        if (bFullW) %-- full W case (not stable)
            wk = W(k,k);
            utB = utB - qold'*W(:,k) + qold(k)*wk;
        else
            wk = W(k);
        end
        
        qnew(k)=softThreshold(utB,lambda1/2)/wk;

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
    if verbose
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

if verbose
    fprintf('#iterations: %d  -- Total Shoots: %d\n',niter-1,(niter-1)*nFea);
end

model.v = theta;
% model.theta_log = theta_log;

%-------------------------------------------------------------------------%

function [scal]=softThreshold(mat,lam)
    %-- compute soft func, both mat and lam are scalars
    scal = sign(mat)*max(abs(mat)-lam, 0);
end

end