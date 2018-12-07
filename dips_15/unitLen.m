function Y = unitLen(X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION
%   Y = unitLen(X)
% Normalize each column (a sample, not feature) of X to unit length
%
% INPUT ARGUMENTS:
%   X[nFea nSmp]: original data
% OUTPUT ARGUMENTS:
%   Y[nFea nSmp]: normalized data 
% (c) 2014 QuangAnh Dang
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[nFea nSmp] = size(X);
lengthX= max(1e-14,full(sum(X.^2)))';
Y= X*spdiags(lengthX.^-.5,0,nSmp,nSmp);

            