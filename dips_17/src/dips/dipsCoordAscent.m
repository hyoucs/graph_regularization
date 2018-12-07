function model = dipsCoordAscent(data, options) % separate lambda1 for parfor

%------------------------------------------------------------------------%
% Minimize the following objective function using coordinate descent
%		R(u) = ||y-Xu||^2 + lambda2*u'Cu + lambda1*|u|
%
% Input:
%     + data: struct data type
%           .X[nSmp nFea]: 	dataset X (or V' in the paper)
%           .y[nSmp 1]: 	mapped point in 1-dim subspace
%           .W[nFea nFea]: 	network adjacency matrix
%     + options: structure
%           .lambda1:   	L1 parameter for sparsness 
%           .lambda2:   	L2 parameter for smoothness
%           .u [nFea 1]: 	initial regression coeff
%           .verbose: 		boolean, used to print out results
%           .nFeaUpd: 		number of features to be randomly updated 
%           (if nFeaUpd < 1, then it is percentage)
%------------------------------------------------------------------------%

	disp('...dipsCoordAscent...');


	% - - - - - - - - - - - - - IMPORT - - - - - - - - - - - - - - - - - 

	% import data
	X = data.X;	% loading features
	y = data.y;	% labels to column vector format
	[nSmp, nFea] = size(X);
	W = full(data.W);	% W maybe sparse		
	W = max(W, W); 	% symmetrize adjacent matrix
	clear data;

	% import parameter for l1 penalty
	lambda1 = 0.1;
	if isfield(options,'lambda1')
	    lambda1 = options.lambda1;
	end

	% import number of features to be updated randomly at each iteration
	nFeaUpd = nFea;
	if isfield(options,'nFeaUpd')
	    if options.nFeaUpd <=1
	        % percentage of no. of features
	        nFeaUpd = floor(options.nFeaUpd*nFea);
	    else
	        % exact number of features
	        nFeaUpd = options.nFeaUpd;
	    end
	end

	% import parameter for quadratic penalty
	lambda2 = 0.1;
	if isfield(options,'lambda2')
	    lambda2 = options.lambda2;
	end

	% import initial point in optimization
	if isfield(options,'u')
	    u = options.u; 
	else
	    u = (X'*X + 0.01*eye(nFea))\(X'*y); % ridge
	    % u = rand(nFea,1); % initialize by random
	end

	% import flag for printing
	verbose = 0;
	if isfield(options,'verbose')
	    verbose = options.verbose; 
	end
	if verbose       %Start the log
	    fprintf('%10s %15s %15s\n','iter','sum(|u|)','sum(|u_{k} - u_{k-1}|)');
	end





	% - - - - - - - - - - - - - DO WORK - - - - - - - - - - - - - - - - - 

	% compute graph Laplacian from adjacency matrix
	D = diag(sum(W,2));
	C = D-W; 
	clear W D;

	% % alternative: symmetric normalized Laplacian
	% D = diag(sum(W,2));
	% C = eye(nFea)-D^(1/2)*W*D^(1/2);
	% clear W D;

	% % alternative: random-walk normalized Laplacian
	% D = diag(sum(W,2));
	% C = eye(nFea)-D^(-1)*W;
	% clear W D;

	% coordinate descent
	maxIter = 500;
	optTol = 1e-5;

	% keep track of coefficient vector u in each iteration
	u_log = zeros(nFea,maxIter);
	u_log(:,1) = u;

	% prepare meta variables
	V = X';
	VVt = V*V';
	Vy = V*y;

	iter 	= 0;
	while iter < maxIter
		u_prev = u;
		% randomly update nFeaUpd features in each iteration 
		rndIdx = randperm(nFea, nFeaUpd);
		for k = 1:nFeaUpd
			t = rndIdx(k);
			a1 = 2 * (VVt(t,t) + lambda2 * C(t,t));
			a2 = 2 * (Vy(t) - VVt(t,:)*u + VVt(t,t)*u(t)) -...
				 2 * lambda2 * (C(t,:) * u - C(t,t) * u(t));
			% update u(t)
			if a2 < - lambda1
				u(t) = (a2 + lambda1) / a1;
			elseif a2 > lambda1
				u(t) = (a2 - lambda1) / a1;
			elseif abs(a2) <= lambda1
				u(t) = 0;
			end	
		end
		% store and print updated coefficient vector
		iter = iter + 1;
		u_log(:,iter+1) = u;
		if verbose
	        fprintf('%10d %15.4e %15.4e\n',iter,sum(abs(u)),sum(abs(u-u_prev)));
	    end
	    % check stopping criterion
	    if sum(abs(u-u_prev)) < optTol
	    	break;
	    end
	end





	% - - - - - - - - - - - - - OUTPUT - - - - - - - - - - - - - - - - - 

	u_log(:,iter:end) = [];

	model.C = C;
	model.nFeaUpd = nFeaUpd;
	model.u = u;
	model.u_log = u_log;

	disp('...Exist dipsCoordAscnet...');

end
