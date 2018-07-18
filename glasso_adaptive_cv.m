function [W,M,M0,rhomax] = glasso_adaptive_cv(Y,rhopath,K,M0)
% K-fold cross-validation (CV) for choosing shrinkage parameter in adaptive graphical
% lasso methods:

%   argmax_M  log(det(M)) - trace(S*M) - rho*weight.*||M||_1

%  using the normal likelihood as objective function for CV.

% Input:  Y: p by n data matrix where dim = p, sample size = n 
%         rhophath: a grid of shrinkage parameter values from which the
%         optimal one is chosen
%         K: number of folds used in the CV
%         M0: input for weight matrix. If not provided, use sample precision matrix
% Output:  
%           W: estimated covariance matrix under optimal shrinkage parameter
%           M: estimated precision matrix under optimal shrinkage parameter
%           M0: weight matrix used
%           rhomax: optimal shrinkage parameter among rhopath
%
%  Ref: Fan Feng Wu 2009 AoAs.
%  Written by Hao Wang & U of South Carolina


[p,n] = size(Y);

if nargin==3
   M0 = inv(Y*Y'/n);
end

   weight = abs(M0).^(-1/2);



k = floor(n/K);

n_train = n -k ;
n_test = k;

loglike = zeros(K,length(rhopath));
for i = 1:K
  disp(['fold ',int2str(i),'/',int2str(K)])
 Y_test = Y(:, (i-1)*k+1:i*k);
 Y_train = Y;
 Y_train( :,(i-1)*k+1:i*k) = [];

 S_train = Y_train*Y_train';
 S_test = Y_test*Y_test';

for j = 1:length(rhopath)
 rho = rhopath(j);   
 [W,M] = glasso_FTH(S_train/n_train,rho.*weight); %L1precisionBCD(S_train/n_train,rho);
 loglike(i,j) = log(det(M)) - trace(S_test/n_test*M); 
end
 
end

[a,b] = max(mean(loglike));

rhomax = rhopath(b);

[W,M] = glasso_FTH(Y*Y'/n,rhomax.*weight);