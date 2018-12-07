function model = dipsLE(data,options,varargin)
%------------------------------------------------------------------------%
%-- Laplacian Embedding, both linear and non-linear projections 
%-- Non-linear: Max y'*L*y s.t. y'*Dp*y=1  
%-- Linear:     v'*X'*L*X*v s.t. v'*X'*Dp*X*v=1  (may add v'*C*v as network topo constraint)
% 
%-- Note: Different sizes of eigVecs y (1*nSmp) and a (1*nFea)
%     Non-linear: entries in y are mapping points
%     Linear: entries in a are linear coeficients 

% To find discriminative subnetworks to classify network samples
% 2 steps:


% seek y s.t. MAX y'(Lb+Aw)y s.t. y'Dwy=1 or equivalently
%               z'(Dw^{-.5}(Lb+Aw)Dw^{-.5})z = z'Lz s.t. z'z=1
%  view L=XX' and decompose X=USV', similar to sparsePCA
%  find sparse a s.t. min ||Xa -y ||^2  + \lambda2 a'Ca + \lambda1 |a|
%  where y' = \sigma*v' by svd

% Input:
%     + data: struct data type
%           .X [nSmp nFea]:     dataset X
%           .gnd [nSmp 1]:      global network states or labels
%           .W [nFea nFea]:     topology sharing among samples
%                 Wij is the frequency of edge ij in all networks
%           .A [nSmp nSmp]:     pairwise similarity (among samples)
%     + options: structure
%           .bLinear: =1 for linear case
%           .k:     number of nearest neighbors, def.=10
%           .beta:  similarity tradeoff (beta*Lb + (1-beta)*Aw), def.=.5
%           .d:     projected subspace dimension, def.= nClass-1
% Output:
%     + model: struct data type
%           .eVecs:             either y or a above. 
%           .eVals:             either y or a above. 
%           .Y [nSmp d]:        samples presented in projected subspace
%           .Dp
%           .L
%------------------------------------------------------------------------%

    %-- beta: affects L1 selection,
    %       better for small beta ~ more emphasis on Aw
    %       if balanced beta, kNN should be large, say 40
    % beta=0.5; %-- Aw is much more important than Lb



% - - - - - - - - - - - - - IMPORT PARAMETERS - - - - - - - - - - - - - - - - - 

    % number of classes
    Class  = unique(data.gnd);
    nClass = length(Class);

    % # features, # samples
    [nSmp, ~] = size(data.X);

    % check existence of parameter options
    if (~exist('options','var'))
        options = [];
    end

    % k: number of nearest neighbors, def.=10
    k = 10;
    if isfield(options,'k')
        k = options.k;
    end

    % bLinear: =1 for linear case
    bLinear = 1;
    if isfield(options,'bLinear')
        bLinear = options.bLinear;
    end

    % beta: similarity tradeoff (beta*Lb + (1-beta)*Aw), def.=.5
    beta = 0.5;
    if isfield(options,'beta')
        beta = options.beta;
    end

    % d: projected subspace dimension, def.= nClass-1
    d=nClass-1;
    if isfield(options,'d')
        d = options.d;
    end




% - - - - - - - - - - - - - DO WORK - - - - - - - - - - - - - - - - - - - - - 

    % - - - - - - - - - STEP 1: COMPUTE SIMILARITY MATRIX (Kp Kn) - - - - - -

    %-- use Euclidean or Pearson for expression profiles
    if ~isfield(options, 'type')
        options.type='Eucl';
    end
    % opt.type='Cosine';
    if isfield(data,'A')
        A = data.A;
        data = rmfield(data,'A');
    else
        A = pwDist(data.X',data.X',options);
    end
    
    % -- generate KNN binary mask
    if options.distOrderAscend % for distance measures, the larger, the less similar
        [~, idx] = sort(A,2,'ascend'); % -- sort each row in a ascend order, 'Eucl'
        idx(:,[1,k+2:end]) = []; % -- remove 1st (selfnode sim) + last largest sim
    else % for similarity measures, the larger, the more similar
        [~, idx] = sort(A,2,'descend'); % -- sort each row in a descend order, 'Cosine'
        idx(:,[1,k+2:end]) = []; % -- remove 1st (selfnode sim) + last smallest sim
    end
    kNN = zeros(size(A));
    for i = 1:nSmp
        kNN(i,idx(i,:)) = 1;
    end

    A = A.*kNN;
    A = max(A,A'); % -- (attention to corr/cosine similarity)


    %-- build 2 meta masks: p(intra-class) and n(inter-class) 
    Ap = zeros(nSmp,nSmp);
    Ap = sparse(Ap);
    for i = 1 : nClass
        classIdx = find(data.gnd==Class(i));
        Ap(classIdx,classIdx) = 1;
    end
    An = ones(nSmp,nSmp)-Ap;

    %-- apply the mask
    Ap = Ap.*A;
    An = An.*A;
    Ap = max(Ap,Ap');
    An = max(An,An');

    % % figure;imshow(full(A),'InitialMagnification',2000);
    % figure; 
    % subplot(1,2,1); spy(Ap);
    % subplot(1,2,2); spy(An);

    % - - - - - - - - - STEP 2: COMPUTE LAPLACIAN (Lp Ln) - - - - - -

    % compute Laplacian Lp of intra-class similarity
    Dp = diag(sum(Ap,2));
    if ~isempty(find(diag(Dp)==0))
        error('Increase kNN to ensure no disconnected same-class nodes!');
    end
    Lp = Dp - Ap; 
    clear Ap; % - - hold Dp

    % compute Laplacian Ln of inter-class similarity
    Dn = diag(sum(An,2));
    Ln = Dn - An; 
    clear An Dn;

    % jointly model inter-class and intra-class similarity
    L = Ln - beta * Lp;
    % L = beta*Ln + (1-beta)*Lp;    % alternative formulation

    % linear case (slower): max a'X(Lb+Aw)X'a s.t. a'XDX'a=1
    if bLinear 
       L = data.X'*L*data.X;
       Dp = data.X'*Dp*data.X;
    end
 


    % - - - - - - - - - STEP 3: SOLVE OPTIMIZATION PROBLEM - - - - - -

    % %-- pseudo inv is worse than using eig
    % [U,S,V]=SVD_CutOff(Dw,.99);
    % invD = V*inv(S)*U';
    % invD = max(invD,invD');
    % [eVecs, eVals] = eigs(invD*L,d); 
    % diag(eVals)
    % Dw*invD

    % %-- better for singular Dw, check algo behind eig and eigs
    % [eVecs, eVals] = eig(L, Dw); % diag(eVals) % work well for non-linear but
    % show for linear case on Brain data


    [eVecs, eVals] = eigs(L, Dp);
    eVals(isnan(eVals)) = 0;
    [eVals, ind] = sort(diag(eVals), 'descend');
    eVecs = eVecs(:,ind(1:min([d size(eVecs, 2)])));
      







% - - - - - - - - - ENCAPSULATION - - - - - - - - - - - - - - - - - -

    model.beta = beta;
    model.k = k;

    model.L = L;
    model.Dp = sparse(Dp);

    model.eVecs = eVecs;
    model.eVals = eVals;

    if bLinear 
        model.Y = data.X * eVecs ;
    else       
        model.Y = eVecs;
    end

    

    % figure;
    %     subplot(1,3,1); PlotX(data.X,data.gnd,'','Original space',''); grid on;
    %     subplot(1,3,2); PlotX(model.Y',data.gnd','','Embedded space 2-dim','');
    %     subplot(1,3,3); stem(abs(model.eVecs(:,1)));

end
