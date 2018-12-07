function 0genData()

%-------------------------------------------------------------------------%
%---------------------  Clemmensen SparseLDA  ----------------------------%
%-------------------------------------------------------------------------%
%-- only 2 Gausses
%%
% Fix stream of random numbers
close all; clear all;clc

s1= RandStream('mrg32k3a','Seed', 50);
s0=RandStream.setGlobalStream(s1);

%   s1 = RandStream.create('mrg32k3a','Seed', 50);
p = 20; % number of variables
nc = 100; % number of observations per class
n = 2*nc; % total number of observations
nGndFea=3;
m1 = 0.6*[ones(nGndFea,1); zeros(p-nGndFea,1)]; % c1 mean
m2 = 0.6*[zeros(nGndFea,1); ones(nGndFea,1); zeros(p-2*nGndFea,1)]; % c2 mean
S = 0.6*ones(p) + 0.4*eye(p); % covariance is 0.6
Yc = [ones(nc,1); 2*ones(nc,1)];

% training data
c1 = mvnrnd(m1,S,nc); % class 1 data
c2 = mvnrnd(m2,S,nc); % class 2 data
X = [c1; c2]; % training data setclc
Y = [[ones(nc,1);zeros(nc,1)] [zeros(nc,1); ones(nc,1)]];

% test data
c1 = mvnrnd(m1,S,nc);
c2 = mvnrnd(m2,S,nc);
X_test = [c1; c2];

figure;PlotX(X(:,1:2)',Yc,'','',''); grid on;

% SLDA parameters
delta = 1e-3;     % l2-norm constraint
stop = -30;       % request 30 non-zero variables
maxiter = 250;    % maximum number of iterations
Q = 1;            % request 1 discriminative direction
convergenceCriterion = 1e-6;

% normalize training and test data
[X mu d] = normalize(X);
X_test = (X_test-ones(n,1)*mu)./sqrt(ones(n,1)*d);

% run SLDA
[B OptScore] = slda(X, Y, delta, stop, Q, maxiter, convergenceCriterion, true);

% Project data onto the sparse directions
DC = X*B;
DC_test = X_test*B;

figure;
subplot(1,2,1);PlotX(DC',Yc,'','',''); grid on;
subplot(1,2,2);PlotX(DC_test',Yc,'','',''); grid on;


% Classification (LDA of projected data)
[class err] = classify(DC, DC, Yc, 'linear');
[class_test] = classify(DC_test, DC, Yc, 'linear');
err_test = sum(Yc ~= class_test)/length(Yc);
fprintf('SLDA: trainErr: %2.1f %%, testErr: %2.1f %%.\n', 100*err, 100*err_test);

[class err] = classify(X, X, Yc, 'linear');
[class_test] = classify(X_test, X, Yc, 'linear');
err_test = sum(Yc ~= class_test)/length(Yc);
fprintf(' LDA: trainErr: %2.1f %%, testErr: %2.1f %%.\n', 100*err, 100*err_test);

% plot sparse discriminative directions for test data
figure;
plot(DC_test(1:nc,1), 1:nc,'ro'), hold on
plot(DC_test((nc+1):2*nc,1), (nc+1):2*nc,'ks')
legend('C_1','C_2','Location','SouthEast')

% Restore random stream
RandStream.setGlobalStream(s0);

feaName = cell(1,p);
for i=1:p
    feaName(i)= {strcat('f',num2str(i))};
end

data.X=X';  data.gnd=Yc';    data.X_test=X_test';    data.feaName = feaName;    data.B=B;    data.OptScore = OptScore;
save('data/2G20Dim3GT.mat','data');


end