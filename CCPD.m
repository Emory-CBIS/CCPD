function[change_points,graph]=CCPD(Y,count)
%%CCPD Toolbox
%Provided by Suprateek Kundu and Jin Ming @ Emory University
%input:Y: Data set, format in T*V*Nsub
%                   T is number of time points, V is number of ROIs, Nsub
%                   is number of subjects
%      count: number of subsampling
%output: change_points: change points for all subjects
%        graph: graph of each bin, format in P*P*Nbin. The max value of
%        Nbin is 20

n=size(Y,1);
p=size(Y,2);
nsub=size(Y,3);
changetime_save=zeros(1,20);
max_jump=20;
jumps = zeros(count,max_jump); % change point matrix
subsam=min(5,round(nsub/10));
nsub1=nsub-subsam;
clsMAT = zeros(count,nsub1); % indices of subjects chosen in a subsample
prob=0.6; % probability for pseudo-probs of subsampling


% subsampling step starts

for rrep = 1:count
    
    cls = datasample(1:nsub,nsub1,'replace',false);
    clsMAT(rrep,:) = cls;
 %   clsMAT_diff(rrep,:) = setdiff(1:nsub,cls);
    Yboot = Y(:,:,cls);

pairwise = zeros(p*(p+1)/2,n);

 for time = 1:n
 tt=0;
 for ii=1:p
     for jj=1:(ii)
     tt = tt+1;    
     pairwise(tt,time) = corr(reshape(Yboot(time,ii,:),[1,nsub1])',reshape(Yboot(time,jj,:),[1,nsub1])') ;
     % check if there is a faster way of computing pairwise corr 
    
     end
 end
 end
 pa = pairwise';
 pa(:,sum(pa,1)==n) = []; %removing the diagonal elements

 
  seg = simpleGFL(pa); % estimating change points based on fused lasso
  jumps(rrep,1:length(seg.jumps)) = seg.jumps;
end


count1 = rrep;
jump_prob = zeros(1,n); % pseudo-probabilities for each time point being a change point

for time=1:n
    for rrep = 1:count1
    jump_prob(time) = jump_prob(time) + mean(ismember(time,jumps(rrep,:)));
    end
    jump_prob(time) = jump_prob(time)/count1;
end


% estimates the changepoints based on pseudo-probs
changetime_est = find(jump_prob(1:(n-1))>prob);
change_points=changetime_est;

changetime_est0 = [0,changetime_est,n];
nbin=length(changetime_est0)-1;
Omega_est = zeros(p,p,length(changetime_est0)-1); % precision matrix for all bins
rho_est = zeros(p,p,length(changetime_est0)-1); % partial corr matrix for all bins

graph=zeros(p,p,nbin);
%% graph estimation
rho1=0.0001:0.01:0.1;
rho2=0.1:0.1:0.8;
rho=[rho1,rho2];
k1=length(rho);
Ome_best=zeros(p,p,nbin);
Ome_all=zeros(p,p,k1,nbin);
bic_all=zeros(k1,nbin);
for kk=1:nbin
Y2=Y((changetime_est0(kk)+1:changetime_est0(kk+1)),:,:);
aa=changetime_est0(kk+1)-changetime_est0(kk);
Y3=zeros(aa*nsub,p);
    for aa2=1:nsub
    Y3((aa2*aa-aa+1):(aa2*aa),:)=Y2(:,:,aa2);
    end
s=cov(Y3);
for k=1:k1
   r1=rho(k);
   [X W opt cputime iter dGap] = QUIC('default', s, r1, 1e-6, 2, 100);
    Ome_all(:,:,k,kk)=X;
    bic_all(k,kk)=-log(det(X))+sum(diag(s*X))+(log(aa*nsub)/(aa*nsub))*sum(sum(abs(X*(tril(X)*1))>0));
end
[M,I]=min(bic_all(:,kk));
Ome_best(:,:,kk)=Ome_all(:,:,I,kk);

end
graph=Ome_best;

end