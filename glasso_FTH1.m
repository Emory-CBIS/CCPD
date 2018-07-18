function [W,M] = glasso_FTH(S,lambda,Max1,Max2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graphical lasso algorithm 
% argmax_M  log(det(M)) - trace(S*M) - lambda* ||M||_1
% Input: S = Y'Y  p*p sample product matrix 
%        lambda>0  shrinkage parameter 
% Output: M estimated precision matrix
%         W = inv(M); estimted covariance matrix
% Ref: Friedman et al 2007 biostatistics

% N.B.  This algorithm may collapse in certain cases, especially when p large, as no positive
% definite contraints are imposed in the glasso algorithm.


% Written by Hao Wang @ U of South Carolina

p = size(S,1);

if isscalar(lambda) == 1
    lambda = lambda*ones(p);
end


idx = find(~eye(p));

avgSoff = mean(S(idx)); % average of off-diagonal elements of empirical covariance matrix S

t = 0.001; % a fixed threshold ; 

%% Initial value
W = S + diag(diag(lambda));

% Maximum number of iterations
Max1 = 30; % across column
Max2 = 30; % Within column, gradient descend

%%
optTol = 1e-4;

for iter1 = 1:Max1
  
  W_old = W;
  
  
  for i = 1:p
      
      
      if i==1  
       ind_noi = [2:p]'; 
      elseif i==p
       ind_noi = [1:p-1]'; 
      else
       ind_noi = [1:i-1,i+1:p]';
      end
      
       V = W(ind_noi,ind_noi);
       

       
       s12 = S(:,i); s12(i) = [];
       w12 = W(ind_noi,i); 
       
      beta = V\w12; 
      
     %% below Pathwise coordinate descent  
      for iter2 = 1:Max2 
       
       beta_old = beta;
       
       for j = 1:p-1           
        
           
        x = s12(j) - V(j,:)*beta + V(j,j)*beta(j);

        
        if x==0
            signx = rand(1)
        else
            signx = sign(x);
        end
        
        beta(j) = max(0,abs(x)-lambda(i,ind_noi(j)))* signx/V(j,j);
        
%         if (W(i,i)-beta'*V*beta<=0)        
%         error('Matrix is not positive definite');
%         end
                
       end
       
       if max(abs(beta_old-beta))< optTol
       break;
       end 
       
      end
      
      w12 = V*beta;
      
      W(i,ind_noi) = w12';
      W(ind_noi,i) = w12;
       
  end
    
      chg = abs(W_old-W);  
      if mean(chg(:)) < t*avgSoff 
       break;
      end 
      
end


M = inv(W);

