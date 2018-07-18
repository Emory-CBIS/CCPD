% function for multi subject fMRI analysis
% nt is number of changepoints
% D is a p X p matrix in the inverse wishart prior. In our method D is
% assumed diagonal D dI_p
% Y is the T X p data matrix
% MAT is the matrix of changetimes given nt changepoints
% numt is the number of timepoints in each bin
function[BIC, numt] = BICcomp(nt,n,nsub,Y,MAT,lambdapath)
cc = size(MAT,1); 
p = size(Y,2);
bins = [zeros(cc,1),MAT,n*ones(cc,1)];
BIC = zeros(cc,nt+1);
numt = zeros(cc,nt+1);
for ii=1:(nt+1)    
    tt1 = bins(:,ii)+ones(cc,1); tt2 = bins(:,ii+1);
    parfor jj=1:cc
        numt(jj,ii) = tt2(jj) - tt1(jj) + 1;
        SS = zeros(p,p);
        for nn=1:nsub
        YY = reshape(Y(tt1(jj):tt2(jj),:,nn),[tt2(jj)-tt1(jj)+1,p]);
        SS = SS + (YY- repmat(mean(YY),[tt2(jj)-tt1(jj)+1,1]) )'*(YY - repmat(mean(YY),[tt2(jj)-tt1(jj)+1,1]) );
        s=cov(YY);
        end
        % [ee, ev] = eig(SS);
        BIC_dum = zeros(length(lambdapath),1); 
     
        for ll=1:length(lambdapath)
        %Omega  = glasso_FTH(SS,lambdapath(ll));
           [X W opt cputime iter dGap] = QUIC('default', s, lambdapath(ll), 1e-6, 2, 100);
           Omega = X;
        BIC_dum(ll) = sum(diag(Omega*SS)) - nsub*numt(jj,ii)*log(det(Omega)) + log(nsub*numt(jj,ii))*(sum(sum(abs(triu(Omega))>0.005))-p);
        end
        BIC(jj,ii) = min(BIC_dum);
    end
    
end
