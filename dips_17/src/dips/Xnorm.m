function [Y]= Xnorm(X,opt)
%------------------------------------------------------------------------
% Normalize each feature (row) of X to have 0-mean and (optional) 1-std.dev. 
%
% INPUT ARGUMENTS:
%   X: DxN matrix of original data in columns
%   opt: not empty then 1-std dev (optional)
% OUTPUT ARGUMENTS:
%   Y: DxN normalized data
%------------------------------------------------------------------------

[D,N]=size(X); 

Y = X-repmat(mean(X')',1,N); %-- 0-mean (repmat: create N mean-cols)
%-- or: Y=bsxfun(@minus,X,mean(X')');
% Y = Y/max(max(abs(Y)));      %-- required for spec. clustering!

%-- one-stand. dev.
if nargin>1
    Y = Y./repmat(std(Y')',1,N); 
end

