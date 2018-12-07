function [beta ] = eNet(data,opt)
% Elasticnet with Shooting lasso optimization
% Input:
%     + data: struct data type
%           .X [nFea nSmp]: normalized 0-mean dataset X
%           .y [1 nSmp]:    normalized 0-mean groundtruth (either continous regression or 
%           .W [nFea nFea]: network topology among nFea
%     + opt: structure
%           .lambda1: L1 regularization
%           .lambda2: L2 regularization
% Output:
%     + model: struct data type
%           .eVecs: either y or a above. 
%           .eVals: either y or a above. 
%           .Y: samples presented in projected subspace (=y in non-linear,
%           = a'X in linear)
%           .Lb = Db-Ab
%           .Aw 
%           .Dw 
% (c) 2014 Quang-Anh Dang - UCSB
%------------------------------------------------------------------------%


X = data.X';
y = data.y';
lambda1 = opt.lambda1;
lambda2 = opt.lambda2;


[nSmp nFea] = size(X);

% network provided
if isfield(data,'W') 
    W = max(data.W,data.W');
else
    W = eye(nFea);
end

%-- construct graph-constrainted Laplacian C=Lc=Dc-Wc
Wc = max(W,W'); 
Dc = full(sum(Wc,2));
C = -Wc; clear Wc
for j = 1 : size(C,1)
	C(j,j) = Dc(j) + C(j,j);
end

[U, D] = ldl(C);
D(D<0)=0;
S = U*D^(.5);

[nSmp nFea] = size(X);
coeff = sqrt(1+lambda2);
Xstar = 1/coeff*[X; sqrt(lambda2)*S']; 
ystar = [y; zeros(nFea,1)];

%-- now call Lasso algorithm
gamma = lambda1/coeff;

% betaStar = LassoShooting(Xstar, ystar, gamma);
betaStar = ShootingAlgo(Xstar, ystar, gamma);

beta = coeff*betaStar; % corrected enet
%beta = betaStar/coeff; % naive enet

end