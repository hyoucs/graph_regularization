function model = regCoordAscent(data,options,lambda1)
% function model = regCoordAscent(data,options,lambda1) % separate lambda1 for parfor
%------------------------------------------------------------------------%
% Max.obj.: ||y-Xa||^2 + lambda2*a'Ca + lambda1*|a|, using coordinate ascent
%
% Input:
%     + data: struct data type
%           .X[nFea nSmp]: dataset X (no need?)
%           .y[1 nSmp]: mapped point in 1-dim subspace
%           .W[nFea nFea]: network topo
%     + options: structure
%           .lambda1:   L1 tradeoff for sparsness 
%           .lambda2:   L2 tradeoff for smoothness
%           .a [nFea 1]: initial regression coeff
%           .verbose: boolean, used to print out results
%           .nFeaUpd: number of features to be randomly updated (nFeaUpd<1,
%           then, it is percentage)
% Output:
%     + model: struct data type
%           .a: sparse vector of coefficients
%           .a_log: log of a till convergence
% Example: (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%

tic
data = c2s(data);
X = data.X;
y = data.y'; % to col vector format
% clear data;
[nFea,nSmp] = size(X);

if isfield(data,'W'),
    Wc = max(data.W,data.W');
else
    Wc = eye(nFea);
end


if (~exist('options','var'))
    options = [];
end

lambda1 = 0.1;
if isfield(options,'lambda1'),
    lambda1 = options.lambda1;
end


% number of features to be updated randomly at each iteration
nFeaUpd = nFea;
if isfield(options,'nFeaUpd'),
    if options.nFeaUpd <=1
        % percentage of no. of features
        nFeaUpd = floor(options.nFeaUpd*nFea);
    else
        % exact number of features
        nFeaUpd = options.nFeaUpd;
    end
end



lambda2 = 0.1;
if isfield(options,'lambda2'),
    lambda2 = options.lambda2;
end

if isfield(options,'a'),
    a = options.a; 
else
    a = (X*X' + 0.01*eye(nFea))\(X*y); % ridge
    % a = rand(nFea,1); % initialize by random
end


verbose = 0;
% if isfield(options,'verbose'),
%     verbose = options.verbose; 
% end

if verbose       %Start the log
    fprintf('%10s %15s %15s\n','iter','sum(|a|)','sum(|a-a_Old|)');
end


%-- construct graph-constrainted Laplacian C=Lc=Dc-Wc
Dc = full(sum(Wc,2));
C = -Wc; clear Wc
for j = 1 : size(C,1)
	C(j,j) = Dc(j) + C(j,j);
end

% %--- not much diff on syn data
% C = diag(1./Dc)*C;%-- random walk Lalp
% % C = diag(1./sqrt(Dc))*C*diag(1./sqrt(Dc));%-- symm Lalp
% C = max(C,C');
% %---


maxIter = 500;
optTol = 1e-5;
a_log = zeros(nFea,maxIter);
a_log(:,1) = a;

% tic;
XXt = X*X'; 
v = diag(XXt); %-- cov. vector 
Xy = X*y;           
iter = 0;
while iter < maxIter
    a_old = a;
%   for j = 1:nFea

    % randomly update nFeaUpd features at each iteration only 
    rdIdx = randperm(nFea);
    for k = 1:nFeaUpd  
        j = rdIdx(k);
        term1 = 2*[v(j) + lambda2*C(j,j)];        
        term2 = 2*[Xy(j) -  XXt(j,:)*a + a(j)*v(j)]...
            - 2*lambda2*[C(j,:)*a - C(j,j)*a(j)]; % XXt(j,:)= X(j,:)*X'


%         term1 = 2*[(1-lambda2)*v(j) + lambda2*C(j,j)];        
%         term2 = 2*(1-lambda2)*[Xy(j) -  XXt(j,:)*a + a(j)*v(j)]...
%             - 2*lambda2*[C(j,:)*a - C(j,j)*a(j)]; 
        
        
        % Update a
        if term2 < -lambda1
            a(j) = (term2 + lambda1)/term1;
        elseif term2 > lambda1
            a(j) = (term2 - lambda1)/term1;
        elseif abs(term2) <= lambda1
            a(j) = 0;
        end
        
        
    end
    
    iter = iter + 1;
    a_log(:,iter+1) = a;

    if verbose
        fprintf('%10d %15.4e %15.4e\n',iter,sum(abs(a)),sum(abs(a-a_old)));
    end
    
    % Check termination
    if sum(abs(a-a_old)) < optTol
        break;
    end
end

a_log(:,iter:end) = [];

model.C = C;
model.nFeaUpd = nFeaUpd;
model.a = a;
model.a_log = a_log;

% fprintf('Elapsed time: %8.2f \n',toc);

end