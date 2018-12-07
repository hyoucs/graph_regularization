function A = vec2mat(v, n)

%------------------------------------------------------------------------%
% For a vector v, reconstruct an n-by-n symmetric matrix A that v 
% contains all entries in A. For example, if v = 
% 	     1
% 	     2
% 	     4
% 	     3
% 	     5
% 	     6
% and n = 4, then
% 	A =
% 	     1     1     2     3
% 	     1     1     4     5
% 	     2     4     1     6
% 	     3     5     6     1
%------------------------------------------------------------------------%


	A = zeros(n,n);
	% map to lower triangular matrix
	A(tril(true(size(A)),-1)) = v;
	% add one to diagonal and symmetrize
	A = A + A' + eye(n);

end