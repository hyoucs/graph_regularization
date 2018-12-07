function v = mat2vec(A)

%------------------------------------------------------------------------%
% For a symmetric matrix A, extract the upper triangle entries 
% into a column vector. For example: 
% 	A =
% 	     1     1     2     3
% 	     1     1     4     5
% 	     2     4     1     6
% 	     3     5     6     1
% Then the vector will be 
% 	v = 
% 	     1
% 	     2
% 	     4
% 	     3
% 	     5
% 	     6
%------------------------------------------------------------------------%

	v = A(triu(true(size(A)),1));

end