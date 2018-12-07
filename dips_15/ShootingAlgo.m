function [w,wp,m] = ShootingAlgo(X,y, lambda,varargin)

%-- DQA: elasticNet also call this function!

[maxIter,verbose,optTol,zeroThreshold] = process_options(varargin,'maxIter',10000,'verbose',0,'optTol',1e-5,'zeroThreshold',1e-4);
[n p] = size(X);

beta = (X'*X + lambda*eye(p))\(X'*y); %-- init beta as OLS


if verbose==2
    fprintf('%10s %10s %15s %15s %15s\n','iter','shoots','abs(bi)','abs(bi-bj)','f(w)');
    wp = beta;
end


A = 2*X'*X;  %-- ~cov matrix
a = diag(A); %-- variance vector
c = 2*X'*y;

m = 0;
while m < maxIter
    beta_old = beta;
    for j = 1:p
        
        aj = a(j);
        cj = c(j) - A(j,:)*beta + beta(j)*aj;
        
        if cj > lambda
            beta(j) = (cj-lambda)/aj;
        elseif cj < -lambda
            beta(j) = (cj+lambda)/aj;
        elseif abs(cj) <= lambda
            beta(j) = 0;
        end
        
    end
    
    m = m + 1;
    
    if verbose==2
        fprintf('%10d %10d %15.4e %15.4e %15.4e\n',m,m*p,sum(abs(beta)),sum(abs(beta-beta_old)),...
            sum((X*beta-y).^2)+lambda*sum(abs(beta)));
        wp(:,m+1) = beta;
    end
    if sum(abs(beta-beta_old)) < optTol
        break;
    end
end
w = beta;
end
