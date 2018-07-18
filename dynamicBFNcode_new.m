

REP = 50;
count = 50;

n=300; p=20; change = [40, 115, 175]; N = length(change)+1;
nsub =60; % Number of subjects
phT = [ones(1,change(1)),2*ones(1,change(2)-change(1)),3*ones(1,change(3)-change(2)),4*ones(1,n-change(3))]; % true allocation indices

TP = zeros(REP,n); TN = zeros(REP,n); FP = zeros(REP,n); FN = zeros(REP,n);
spec = zeros(REP,1); sens = zeros(REP,1);
changetime_save = zeros(REP,10);


for lp = 1: REP
tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Graph generation
%%% Diagonally dominant case incorporating covariate info
OmegaT = zeros(p,p,N); SigTrue = zeros(p,p,N);
probMAT = zeros(p,p,N);
common = [0,0];
frac = round((p*(p-1)/2)/5);
cls = unidrnd(2,1);
edge = zeros(p,p);
diffmat2 = zeros(frac,2,N);

for nn=1:N
% generate the data for the first group
if(nn==1)
prob = 0.25;
    for i=1:p
    for j=(i+1):p 
    probMAT(i,j,nn) = prob; probMAT(j,i,nn) = prob;
    % has an edge with probability prob
    edge(i,j) = binornd(1,prob);
    edge(j,i) = edge(i,j);
    OmegaT(i,j,nn) = edge(i,j)*unifrnd(-1,1); OmegaT(j,i,nn) = OmegaT(i,j,nn);
    OmegaT(i,i,nn) = 1;
    
    end
    end

% generate the data for the rest of the groups
else
    % row, col
    [v1,v2]=find(triu(edge) - eye(p) >0); % index of present edges
    V1 = [v1,v2];
    [vo1,vo2] = find( (ones(p,p) - triu(edge))==1);  %index of absent edges
    V0 = [vo1,vo2];
    % random samples of size frac/2 from 1:length(V1 or V2)
    diffrnd1 = randsample(size(V1,1), round(frac/2)); 
    diffrnd0 = randsample(size(V0,1), round(frac/2));
    
    edge1 = edge;
    for ii=1:round(frac/2)
        % CHANGE IS HERE
        index1 = diffrnd1(ii);
        index2 = diffrnd0(ii);
        % set some of the ones with edges to zero
        edge1(v1(index1),v2(index1)) = 0; edge1(v2(index1),v1(index1)) = edge1(v1(index1),v2(index1));
        % set some that were 0 to have edges
        edge1(vo1(index2),vo2(index2)) = 1; edge1(vo2(index2),vo1(index2)) =     edge1(vo1(index2),vo2(index2));
    end
    for i=1:p
    for j=(i+1):p
        % get strength of edges
        OmegaT(i,j,nn) = edge1(i,j)*unifrnd(-1,1); OmegaT(j,i,nn) = OmegaT(i,j,nn);
        % make diags 1
        OmegaT(i,i,nn) = 1;
    end
    end
    edge = edge1;
    end
 end
    
  
for nn=1:N
[ee,ev] = eig(reshape(OmegaT(:,:,nn),[p,p]));
OmegaT(:,:,nn) = OmegaT(:,:,nn) + (0.1 -min(diag(ev)))*eye(p);


SigTrue(:,:,nn) = inv(reshape(OmegaT(:,:,nn),[p,p]));
SigTrue(:,:,nn) = (reshape(SigTrue(:,:,nn),[p,p]) + reshape(SigTrue(:,:,nn),[p,p])')/2;
end

spike = zeros(n,1);
spike(20) = 4; spike(100) = 4; spike(200)=4;
% Training Data Generation
Y = zeros(n,p,nsub);
for ii=1:n
    for jj=1:nsub
    if(ii<=change(1))
        Y(ii,:,jj) =  spike(ii) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,1),[p,p]));
    elseif ((ii>change(1))&&(ii<=change(2)))
        Y(ii,:,jj) =  spike(ii) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,2),[p,p]));
    elseif ((ii>change(2))&&(ii<=change(3)))
        Y(ii,:,jj) =  spike(ii) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,3),[p,p]));
    elseif (ii>change(3))
        Y(ii,:,jj) =  spike(ii) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,4),[p,p]));
    end
    end
   
end


% spikes: 80, 150, 300 with real change points:40, 115, 175
% Y(80,:,:)=Y(80,:,:)+4;
% Y(150,:,:)=Y(150,:,:)+4;
% Y(300,:,:)=Y(300,:,:)+4;

Yavg = zeros(n,p); % averaged values over number of subjects
for ii=1:n
    for jj=1:p
        Yavg(ii,jj) = mean(Y(ii,jj,:));
    end
end



% Test Data Generation
nsub_test = 30;
Ytest = zeros(n,p,nsub_test);
for ii=1:n
    for jj=1:nsub_test
    if(ii<=change(1))
        Ytest(ii,:,jj) =  mvnrnd(zeros(1,p),reshape(SigTrue(:,:,1),[p,p]));
    elseif ((ii>change(1))&&(ii<=change(2)))
        Ytest(ii,:,jj) = -5*ones(1,p) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,2),[p,p]));
    elseif ((ii>change(2))&&(ii<=change(3)))
        Ytest(ii,:,jj) = 0*ones(1,p) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,3),[p,p]));
    elseif (ii>change(3))
        Ytest(ii,:,jj) =  5*ones(1,p) + mvnrnd(zeros(1,p),reshape(SigTrue(:,:,4),[p,p]));
    end
    end
   
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% screening step, gives us likely changepoints after screening out
% unnecessary ones
% [sall]=screen1(Yavg,10,1);
% 
% pairwise = zeros(p*(p+1)/2,n);
% for time = 1:n
% tt=0;
% for ii=1:p
%     for jj=1:(ii)
%     tt = tt+1;    
%     pairwise(tt,time) = corr(reshape(Y(time,ii,:),[1,nsub])',reshape(Y(time,jj,:),[1,nsub])') ;
%     end
% end
% end
% 
% pa = pairwise';
% screen2 = screen1(pa,10,3.9);
% 
% % possible changepoints, NN is the number of bins
numchange = [2,3,4,5]; NN = max(numchange)+1;

max_jump = 10;
jumps = zeros(count,max_jump);
for rrep = 1:count
    
    nsub1 = nsub-5;
    cls = datasample(1:nsub,nsub1,'replace',false);
    Yboot = Y(:,:,cls);

pairwise = zeros(p*(p+1)/2,n);
 for time = 1:n
 tt=0;
 for ii=1:p
     for jj=1:(ii)
     tt = tt+1;    
     pairwise(tt,time) = corr(reshape(Yboot(time,ii,:),[1,nsub1])',reshape(Yboot(time,jj,:),[1,nsub1])') ;
     end
 end
 end
 
 pa = pairwise';
 pa(:,sum(pa,1)==n) = []; %removing the diagonal elements
% s2creen = screen1(0.5*log((1+pa)./(1-pa)),10,3.5); % screening timepoints based on Fisher transform for
 %the pairwise correlations
 
  seg = simpleGFL(pa);
  [lp;rrep;seg.jumps]
  jumps(rrep,1:length(seg.jumps)) = seg.jumps;
  
end

jump_prob = zeros(1,n);
for time=1:n
    for rrep = 1:count
    jump_prob(time) = jump_prob(time) + mean(ismember(time,jumps(rrep,:)));
    end
    jump_prob(time) = jump_prob(time)/count;
end

Ymean = mean(Y,3);
changetime_est = find(jump_prob(1:(n-1))>0.8);
changetime_save(lp,1:length(changetime_est)) = changetime_est;
changetime_est0 = [0,changetime_est,n];
Omega_est = zeros(p,p,length(changetime_est0)-1);
rho_est = zeros(p,p,length(changetime_est0)-1);

ph_est = zeros(n,1);
for ii=2:length(changetime_est0)
    ph_est(changetime_est0(ii-1)+1:changetime_est0(ii)) = ii-1;
    Yvec = [];
    for sub = 1:nsub
    Yvec = [Yvec;Y(changetime_est0(ii-1)+1:changetime_est0(ii),:,sub) - Ymean(changetime_est0(ii-1)+1:changetime_est0(ii),:)];
    end
[Sigma,Omega,rhomax] = glasso_cv(Yvec',linspace(0,1,100),5);
Omega_est(:,:,ii-1) = Omega;
rho_est(:,:,ii-1) = 2*eye(p) - Omega_est(:,:,ii-1)./sqrt(diag(Omega_est(:,:,ii-1))*diag(Omega_est(:,:,ii-1))');
end

 for time = 1:n
 TN(lp,time)=0; TP(lp,time) = 0; FP(lp,time) = 0; FN(lp,time) = 0; 
 for j=1:p
     ind = 1:p; ind(j) = [];
     TN(lp,time) = TN(lp,time) + sum( (abs(rho_est(j,ind,ph_est(time)))<=0.01).*(abs(OmegaT(j,ind,phT(time))) ==0) );
     TP(lp,time) = TP(lp,time) + sum( (abs(rho_est(j,ind,ph_est(time)))>0.01).*(abs(OmegaT(j,ind,phT(time))) >0) );
     FP(lp,time) = FP(lp,time) + sum( (abs(rho_est(j,ind,ph_est(time)))>0.01).*(abs(OmegaT(j,ind,phT(time))) ==0) );
     FN(lp,time) = FN(lp,time) + sum( (abs(rho_est(j,ind,ph_est(time)))<=0.01).*(abs(OmegaT(j,ind,phT(time))) >0) );
 end
 end
 
 spec(lp) = sum(TN(lp,:))/sum(TN(lp,:)+FP(lp,:)); 
 sens(lp) = sum(TP(lp,:))/sum(TP(lp,:)+FN(lp,:));
 [lp, spec(lp), sens(lp)]
 toc
end
