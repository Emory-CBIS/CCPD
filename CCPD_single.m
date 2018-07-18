%CCPD_single
%Provided by Suprateek Kundu and Jin Ming @ Emory University
%input:Y: Data set, format in T*V*Nsub
%                   T is number of time points, V is number of ROIs, Nsub
%                   is number of subjects
%output: change_points: change points for each single subject
%        graph: graph of each bin, format in P*P*Nbin*Nsub. The max value of
%        Nbin is 20


function[change_points,graph]=CCPD_single(Y)
n=size(Y,1);
p=size(Y,2);
nsub=size(Y,3);
count=nsub*(nsub-1)/2;
max_jump=20;
jumps = zeros(count,max_jump); % change point matrix


% subsampling step starts
clsMAT=nchoosek(1:nsub,2);
for rrep = 1:count
    cls = clsMAT(rrep,:);
 %   clsMAT_diff(rrep,:) = setdiff(1:nsub,cls);
    Yboot = Y(:,:,cls);

pairwise = zeros(p*(p+1)/2,n);
nsub1=2;
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

%% change point estimation for each single subject

sub_jum3_all=zeros(nsub,nsub*20);  %all possible change points
jum_bic2_all=zeros(2,20,nsub);  %BIC bootstrap result
jum_pop3_all=zeros(2,20,nsub);  %grouping results without BIC
dist_max=4;
for sub=1:nsub
index=zeros(count,1);
for i=1:count
index(i,1)=ismember(sub,clsMAT(i,:));
end
% tbl_jum2: all possible points for sub
sub_jum=jumps(index==1,:);
sub_jum2=nonzeros(sub_jum);
sub_jum3=sort(sub_jum2)';
sub_jum3=sub_jum3(sub_jum3<=n-20);  %delete extreme value of change points
sub_jum3=sub_jum3(sub_jum3>=20);

%graph create
sub_jum3_all(sub,1:size(sub_jum3,2))=sub_jum3;

tbl_jum=tabulate(sub_jum3);
tbl_jum2=tbl_jum(tbl_jum(:,2)~=0,:); 

%delete those show up 1 times change happend here!!!!!!
tbl_jum2=tbl_jum2(tbl_jum2(:,2)~=1,:);

% jum_pop2/3 grouping results
jum_pop=zeros(size(tbl_jum2,1),size(tbl_jum2,1)); %the grouping of all possible change points
jum_pop3=zeros(size(tbl_jum2,1),2); %the value of each group and sum of all prob
jum_prob=zeros(1,size(tbl_jum2,1));
if tbl_jum2(1,1) ~= 300
for i=1:size(tbl_jum2,1)
   if size(tbl_jum2,1)==1
            jum_pop(1,1)=tbl_jum2(1,1);
            jum_pop3(1,1)=tbl_jum2(i,1);
            jum_pop3(1,2)=tbl_jum2(i,2);
   else 
   if i==1 
       n_row=1;
       n_col=1;
       jum_pop(1,1)=tbl_jum2(1,1);
       jum_prob(1,1)=tbl_jum2(1,2);
   elseif i<size(tbl_jum2,1)
       if tbl_jum2(i,1)-tbl_jum2(i-1,1)<=dist_max
           n_col=n_col+1;
           jum_pop(n_row,n_col)=tbl_jum2(i,1);
           jum_prob(1,n_col)=tbl_jum2(i,2);
       else
           [~,ind]=max(jum_prob);
           jum_pop3(n_row,1)=jum_pop(n_row,ind);
           %[ind,~]=max(jum_prob);
           %jum_pop3(n_row,1)=round(mean(jum_pop(n_row,jum_prob(1,:)==ind)));
           %jum_pop3(n_row,1)=round(jum_pop(n_row,:)*jum_prob'/sum(jum_prob));
           jum_pop3(n_row,2)=sum(jum_prob);
           jum_prob=zeros(1,size(tbl_jum2,1));
           n_col=1;
           n_row=n_row+1;
           jum_pop(n_row,n_col)=tbl_jum2(i,1);
           jum_prob(1,n_col)=tbl_jum2(i,2);
       end
   else
       if tbl_jum2(i,1)-tbl_jum2(i-1,1)<=dist_max
           n_col=n_col+1;
           jum_pop(n_row,n_col)=tbl_jum2(i,1);
           jum_prob(1,n_col)=tbl_jum2(i,2);
           [~,ind]=max(jum_prob);
           jum_pop3(n_row,1)=jum_pop(n_row,ind);
           jum_pop3(n_row,2)=sum(jum_prob);
       else
           [~,ind]=max(jum_prob);
           jum_pop3(n_row,1)=jum_pop(n_row,ind);
           n_col=1;
           n_row=n_row+1;
           jum_pop3(n_row,1)=tbl_jum2(i,1);
           jum_pop3(n_row,2)=tbl_jum2(i,2);
           jum_pop(n_row,1)=tbl_jum2(i,1);
       end
   end
end
end

jum_pop2 = jum_pop(all(jum_pop==0,2)==0,:);
n_col=size(tbl_jum2,1)-min(sum(jum_pop2==0,2));
jum_pop2=jum_pop2(:,1:n_col);
for i=1:size(jum_pop2,1)
jum_pop3(i,2)=sum(tbl_jum2(ismember(tbl_jum2(:,1),jum_pop2(i,:)),2));
end
jum_pop3= jum_pop3(all(jum_pop3==0,2)==0,:);

jum_pop3_all(:,1:size(jum_pop3,1),sub)=jum_pop3';
%BIC method to re-decrease the number of FP
jum_bic = jum_pop3; %BIC reduction result
n_bic = size(jum_bic,1);
[~,bic_ind]=sort(jum_bic(:,2),'descend');
jum_bic_test = jum_bic(bic_ind,:);
jum_bic2=[];
n_perm=100;
if n_bic <=3
    jum_bic2=jum_bic;
else 
     jum_bic2=jum_bic_test(jum_bic_test(:,2)>jum_bic_test(4,2),:);
    jum_bic_test = jum_bic_test(jum_bic_test(:,2)<=jum_bic_test(4,2),:);
    n_test = size(jum_bic_test,1);
    bic_test=zeros(1,n_test);
    for i=1:n_test
     changetime_est0 = sort([0,jum_bic2(:,1)',jum_bic_test(i,1),n]);
       ll = length(changetime_est0)-1;
       Yboot =   Y(:,:,sub) ;
       bic_test(i) = sum(BICcomp2(ll-1,n,1,Yboot,changetime_est0(2:ll),linspace(0,0.1,10)) );
    end
    [~,bic_ind2]=sort(bic_test);
    jum_bic_test2=jum_bic_test(bic_ind2,:);
    for i=1:n_test
        pos_ind=sum(jum_bic2(:,1)<jum_bic_test2(i,1));
        changetime_est0 = sort([0,jum_bic2(:,1)',jum_bic_test(i,1),n]);
        nt=changetime_est0(pos_ind+3)-changetime_est0(pos_ind+1);
        Yboot=Y(changetime_est0(pos_ind+1)+1:changetime_est0(pos_ind+3),:,sub);
        changetime_est1 = changetime_est0((pos_ind+1):(pos_ind+3))-changetime_est0(pos_ind+1);
        perm_test=sum(BICcomp2(1,nt,1,Yboot,changetime_est1(2),linspace(0,0.1,10)) );
            bic_perm=zeros(1,n_perm);
      for j=1:n_perm
            Yboot2=Yboot(randperm(nt),:);
            bic_perm(j) = sum(BICcomp2(1,nt,1,Yboot2,changetime_est1(2),linspace(0,0.1,10)) );
      end
        if perm_test < quantile(bic_perm,0.2)
        jum_bic2(size(jum_bic2,1)+1,:)=jum_bic_test2(i,:);
        [~,bic_ind3]=sort(jum_bic2(:,1));
        jum_bic2=jum_bic2(bic_ind3,:);
        end
    end
end
jum_bic2_all(:,1:size(jum_bic2,1),sub)=jum_bic2';
end
end
change_points=reshape(jum_bic2_all(1,:,:),[nsub 20])';

%% graph estimation
nbin_max=20;  %max number of bins for each subject

graph=zeros(p,p,nbin_max,nsub);
rho1=0.0001:0.01:0.1;
rho2=0.1:0.1:0.8;
rho=[rho1,rho2];
k1=length(rho);
Ome_best=zeros(p,p,nbin_max,nsub);
Ome_all=zeros(p,p,k1,nbin_max,nsub);
bic_all=zeros(k1,nbin_max,nsub);
for jj=1:nsub
    cpt=jum_bic2_all(1,:,jj);
    cpt2=nonzeros(cpt)';
    cpt2=sort(cpt2);
    bin=[0,cpt2,n];
    nbin=length(bin)-1;
    Y2=Y(:,:,jj);
    changetime_est0=bin;
for kk=1:nbin
Y3=Y2((changetime_est0(kk)+1:changetime_est0(kk+1)),:);
aa=changetime_est0(kk+1)-changetime_est0(kk);
s=cov(Y3);
for k=1:k1
   r1=rho(k);
   [X W opt cputime iter dGap] = QUIC('default', s, r1, 1e-6, 2, 100);
    Ome_all(:,:,k,kk,jj)=X;
    bic_all(k,kk,jj)=-log(det(X))+sum(diag(s*X))+(log(aa*nsub)/(aa*nsub))*sum(sum(abs(X*(tril(X)*1))>0));
end
[M,I]=min(bic_all(:,kk,jj));
Ome_best(:,:,kk,jj)=Ome_all(:,:,I,kk,jj);

end
end
graph=Ome_best;

end
