function beta = elasticNet(X, y, lambda1, lambda2, lassoSolver, varargin)
% elastic net is L1 + L2 regularized least squares regression
% We assume X is standardized and y is centered
% lassoSolver must have this interface:
 %  w = lassoSolver(X, y, lambda1, varargin{:})
% See also elasticNetPath 

% This file is from pmtk3.googlecode.com


if nargin < 5 %-- if lassoAlg is not provided, choose ShootingAlgo as default
  lassoSolver = @LassoShooting;
end

% %-- test with Lasso and Ridge
lambda2=0; %-- L1 only
% % lambda1=0; %-- L2 only
% %-- end test




[n d] = size(X);

%-- DQA: 3 steps below is to change dataset X = (X; S*I), Y =(Y; zero) as by eNet
S = sqrt(1+lambda2);
Xstar = 1/S*[X; sqrt(lambda2)*eye(d,d)]; 
ystar = [y; zeros(d,1)];

%-- now call Lasso algorithm
gamma = lambda1/S;
betaStar = lassoSolver(Xstar, ystar, gamma, varargin{:});
beta = S*betaStar; % corrected enet
%beta = betaStar/S; % naive enet

end
