function [PWdist] = pwDist(X,Y,opt)
%-------------------------------------------------------------------------%
% Compute Gauss/Eucl/Cosine/Pearson Correlation btw points of Y to ALL points of X
% To avoid "out of mem": Y should contain 1 point/vector
% 
% Note:
%     + (comment on values of each type here)
%     + For text data: cosine should be used as by HanJaweiBook,pg397- Vector
%       Objects, it is a nonmetric similarity. Since doc-vectors are non-negative, 
%       cos() is btw 0 (independent) and 1 (most similar). 
%     + For other cases: vectors are both pos/neg entries, not good as
%       cos() is btw -1 and 1. As cos(a,b) = 0 (a,b are independent) is even 
%       more similar than cos(a,b) < 0!
% 
%     + Gene expression profiles: Euclidean distance and Pearson
%       correlation are most commonly used
% Input:
%     + X: set of points in nFea*nSmp1 matrix
%     + Y: set of points in nFea*nSmp2 matrix 
%     + opt:
%         - .type: Eucl (default)/Gauss/Cosine/Pearson
%         - .sigma : for Gauss distance
%               = 0, computing t from X  (see constructW.m from DC).  
% Output: 
%   similarity matrix size of nSmp2*nSmp1 (btw any yj in Y to xi in X)
%     + PWdist: pairwise distance 
% Note: squared E-dist does not affect kNN searching!
%       E(i,j) is btw [0,+inf) AND  G(i,j) is btw [1,0] 
%       G(i,j) -->0 as E(i,j)--> +inf
%       (i,i): E=0 is min while G=1 (or non-exp term if full G used) is max
% Example: Compare:
%        X= rand(10,3); %-- 3 points in R^10
%           opt=[];
%           opt.type='Gauss';  [PWdist] = PWdistance(X,X,opt)   
%           opt.type='Gauss';  opt.sigma=1; [PWdist] = PWdistance(X,X,opt)   
%           opt.type='Eucl';    [PWdist] = PWdistance(X,X,opt)   
%           opt.type='Cosine';  [PWdist] = PWdistance(X,X,opt)   
%           opt.type='Pearson'; [PWdist] = PWdistance(X,X,opt)   
%-------------------------------------------------------------------------%


    [nFea1 nSmp1] = size(X);
    [nFea2 nSmp2] = size(Y);

    if nFea1 ~= nFea2
        error('nFea of X and Y are not matched');
    else
        nFea = nFea1;
        clear nFea1 nFea2
    end

    if ~isfield(opt,'type'),        
        error('A Type of distance (Eucl/Gauss/Cosine/Pearson) should be provided');
    end

    PWdist=zeros(nSmp2,nSmp1); 

    switch opt.type
        case 'Eucl'
            for j=1:nSmp2
                PWdist(j,:) = sqrt(sum((Y(:,j)*ones(1,nSmp1)-X).^2));
                if mod(j,10) == 0
                    fprintf('Pairwise Dist: %d in %d\n',j,nSmp2);
                end
            end
        case 'Gauss'
            E = zeros(nSmp2,nSmp1); 
            for j=1:nSmp2
                E(j,:)=sum((Y(:,j)*ones(1,nSmp1)-X).^2);
            end
            if isfield(opt,'sigma') && (opt.sigma> 0)
                sigma=opt.sigma;
            else %-- compute kernel width from data
                %-- C1
                sigma = mean(diag(cov(X')));
                sigma = 1.06*sigma*(nSmp1)^(-.25); %-- see Vu DASFAA''11, pg10 [17]
                %-- C2 (DC's costructW.m)
            end    
            %-- E=dist.^(0.5);
            PWdist=exp(-E/(2*sigma)); 
        case 'Cosine'
            %-- unitLen ensure each sample has 1-length
            %-- identical to DC if Xnorm is not used
            PWdist =unitLen(Y)'*unitLen(X);
        case 'Pearson'
            for j=1:nSmp2
                PWdist (j,:)=corr(Y(:,j),X);
            end
    end
            
end


function Y = unitLen(X)

    [nFea nSmp] = size(X);
    lengthX= max(1e-14,full(sum(X.^2)))';
    Y= X*spdiags(lengthX.^-.5,0,nSmp,nSmp);

end