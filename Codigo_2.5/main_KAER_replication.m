%%% This code replicates Figure 3 in Kilian (AER 2009), "Not All Oil Price Shocks 
%%% Are Alike: Disentangling Demand and Supply Shocks in the Crude Oil Market"
%%% using the BH algorithm and includes posterior distributions of key
%%% parameters;
%%% See Figures 1-2 in: Baumeister and Hamilton, "Structural Interpretation 
%%% of Vector Autoregressions with Incomplete Identification: Revisiting 
%%% the Role of Oil Supply and Demand Shocks" (December 2017)

%%% Trivariate Oil Market VAR with uninformative priors on lagged coefficients 

%%% Identification of 3 shocks with zero restrictions: 
%%% oil supply shock, aggregate demand shock, oil-specific demand shock

%%% By Christiane Baumeister and James D. Hamilton (Sept 2014)

clear; 
clc;

ndraws=2e5;     %number of MH iterations (paper uses 2e6)
nburn=1e5;      %number of burn-in draws (paper uses 1e6)

% NOTE:
% The parameter xsi is the tuning parameter for the MH step and has to be
% re-set for each application such that it implies a 30-35% acceptance ratio.
% If acceptance ratio is too high, make xsi bigger; 
% if acceptance ratio is too low, make xsi smaller.
xsi=1.2^2;       %tuning parameter to achieve a target acceptance rate of around 35%

select_sample=2; %1: updated sample; 2: KM12 dataset from JEEA; 3: K dataset from AER

nlags=24;        %number of lags
hmax=18;         %impulse response horizon (17 months)
ndet=1;          %number of deterministic variables 
                 %(1: constant; 2: constant and time trend)

if select_sample==1
    load data_k3.txt                %sample: 1973M1 to 2014M1
    % column 1: global crude oil production (in million barrels per day)
    % column 2: global real economic activity index (Kilian, AER 2009)
    % column 3: real refiners' acquisition cost of crude oil imports (deflated with US CPI)

    % transformations of variables
    qo=lagn(100*log(data_k3(:,1)),1);
    rea=data_k3(2:end,2);
    rpo=100*log(data_k3(2:end,3));
    yy=[qo rea rpo];
    yy=yy(1:428,:);

    time=(1973+1/12:1/12:2014)';     %sample period: 1973M2 to 2014M1

elseif select_sample==2
       load kiliandata.txt           %sample: 1973M2 to 2008M9 (from JEEA)
       yy=kiliandata;
       time=(1973+1/12:1/12:2008+8/12)';
       
elseif select_sample==3
       load data_AER                 %sample: 1973M2 to 2006M12
       yy=y;
       clear y
       time=(1973+1/12:1/12:2006+11/12)';
end

s=size(time,1);
n=size(yy,2);                 %number of endogenous variables
                              
% Get data matrices
[X,y,T]=getXy(yy,ndet,nlags);
yyy=y;
xxx=X;

seednumber=140778;
rand('seed',seednumber);
randn('seed',seednumber);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%     ALGORITHM FOR GETTING POSTERIORS OF A, B and D    %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1a: Set parameters of the prior distributions for impact coefficients (A) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
bounds = 5;
%x1=-bounds:.01:0;       %grid for negative parameters
%y1=0:.01:bounds;        %grid for positive parameters
z1=-bounds:.01:bounds;   %grid for parameters where no sign is imposed a priori

% eta(1): short-run oil production elasticity of global demand (no sign)
c_eta1 = 0;
sigma_eta1 = 100;   
nu_eta1 = 3;
prior_eta1 = student_prior(z1,c_eta1,sigma_eta1,nu_eta1);

% gamma_world(2): income elasticity of oil demand 
c_gamma2 = 0;
sigma_gamma2 = 100;   
nu_gamma2 = 3;
prior_gamma2 = student_prior(z1,c_gamma2,sigma_gamma2,nu_gamma2);

% beta(d): short-run inverse price elasticity of oil demand 
c_beta = 0;     
sigma_beta = 100;   
nu_beta = 3;
prior_beta = student_prior(z1,c_beta,sigma_beta,nu_beta);

% Set arbitrary initial values for elements of A 
A_old=[c_eta1; c_beta; c_gamma2];   %prior mode of elements in A
c=size(A_old,1);                    %number of parameters to be estimated in A

% Construct full matrix A (equation (43))
% A=[1 0 0; -A_old(1,1) 1 0; -A_old(2:3,1)' 1];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1b: Set informative priors on lagged coefficients (B) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute standard deviation of each series residual via an OLS regression
% to be used in setting the prior (here: AR(4))
[s11,uhat1]=sd_prior(yy(:,1),nlags);
[s22,uhat2]=sd_prior(yy(:,2),nlags);
[s33,uhat3]=sd_prior(yy(:,3),nlags);

% See Doan (2013) for choices of values of hyperparameters (Minnesota prior)
lambda0=10^9;    %overall confidence in prior (smaller lambda0 implies greater weight to RW)
lambda1=1;       %confidence on higher-order lags (lambda1 = 0 gives all lags equal weight)
lambda2=1;       %confidence in other-than-own lags 
lambda3=100;     %tightness of constant term 

% Specify the prior mean of the coefficients of the 3 equations of the VAR
% and their prior covariance

% EXPECTED VALUE OF REDUCED-FORM COEFFICIENTS
eta=[eye(n) zeros(n,(n*nlags+1)-n)];

% PRIOR COVARIANCE 
SS=[s11;s22;s33];
M=getM(lambda0,lambda1,lambda2,lambda3,SS,nlags,n);
M1=M(:,:,1);
M2=M(:,:,2);
M3=M(:,:,3);

% Compute summary statistics of the observed data:
Syy=yyy'*yyy;
Sxx=xxx'*xxx;
Sxy=yyy'*xxx;
omega_hatT=(Syy-Sxy*inv(Sxx)*Sxy')/T;

% Compute M_star(i) 
M_star1=inv(Sxx+inv(M1));
M_star2=inv(Sxx+inv(M2));
M_star3=inv(Sxx+inv(M3));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1c: Set priors for inverse of diagonal elements (D) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The mean is calibrated on diagonal elements in omega
uhat=[uhat1 uhat2 uhat3];
S=uhat'*uhat/T;

kappa=0.5;
kappastar=kappa+(T/2); % posterior kappastar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JIM's PROPOSAL to deal with scaling issue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixed parameter values
param=[c_eta1;sigma_eta1;nu_eta1;c_gamma2;sigma_gamma2;nu_gamma2; ...
       c_beta;sigma_beta;nu_beta;T;vec(omega_hatT);kappastar;kappa];
   
f_anon = @(theta_hat)post_val(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M_star1,M_star2,M_star3,S);

% find posterior mode
theta_zero = A_old;     %mode of prior distributions
options = optimset('LargeScale','off');
[theta_max,val_max,exitm,om,gm,HM] = fminunc(f_anon,theta_zero,options);
PH=chol(inv(HM))';

%start MH algorithm with theta_hat
A_old=theta_max;
A=[1 0 0; -A_old(1,1) 1 0; -A_old(2:3,1)' 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Set the variance of the candidate generating density (P) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
W=xsi*eye(c);     %variance of RW-MH  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Evaluate posterior at starting value for A:  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega=A*S*A';
taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);

% Evaluate prior p(A) at old draw
prior_e = student_prior(A_old(1,1),c_eta1,sigma_eta1,nu_eta1);
prior_c = student_prior(A_old(3,1),c_gamma2,sigma_gamma2,nu_gamma2);
prior_b = student_prior(A_old(2,1),c_beta,sigma_beta,nu_beta);

% Compute posterior value at candidate draw
log_priors=log(prior_e)+log(prior_c)+log(prior_b);
up=log_priors+T/2*log(det(A*omega_hatT*A'));
down=kappastar*log((2/T)*taustar1) ...
    +kappastar*log((2/T)*taustar2) ...
    +kappastar*log((2/T)*taustar3);
posteriorOLD=up-down;

% RW-MH algorithm 
naccept=0;
count=0;

% Store posterior distribution (after burn-in)
A_post=zeros(c,ndraws-nburn);
IRF=zeros(n,n,hmax,ndraws-nburn);

while count<ndraws 
      count=count+1;
      if (count/10000) == floor(count/10000)
          count
      end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % STEP 4a: Generate draw for A from the RW candidate density %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_new=A_old+chol(W)'*PH*randn(c,1)/sqrt(0.5*(randn(1)^2 + randn(1)^2));    % fat tails
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4b: Evaluate posterior at new draw %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
       % Evaluate prior p(A) at new draw
       prior_e = student_prior(A_new(1,1),c_eta1,sigma_eta1,nu_eta1);
       prior_b = student_prior(A_new(2,1),c_beta,sigma_beta,nu_beta);
       prior_c = student_prior(A_new(3,1),c_gamma2,sigma_gamma2,nu_gamma2);
       
       % Construct full matrix A 
       A=[1 0 0; -A_new(1,1) 1 0; -A_new(2:3,1)' 1];
       omega=A*S*A';
       taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
       taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
       taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
           
       % Compute posterior value at new candidate draw
       log_priors=log(prior_e)+log(prior_b)+log(prior_c);
       up=log_priors+T/2*log(det(A*omega_hatT*A'));
       down=kappastar*log((2/T)*taustar1) ...
           +kappastar*log((2/T)*taustar2) ...
           +kappastar*log((2/T)*taustar3);
       posteriorNEW=up-down; 
         
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 5: Compute acceptance probability %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       accept=min([exp(posteriorNEW-posteriorOLD);1]);
       u=rand(1);                      %draw from a uniform distribution
       if u<=accept
          A_old=A_new;                %we retain the new draw
          posteriorOLD=posteriorNEW;
          naccept=naccept+1;          %count the number of acceptances
       end
         
       
    if count>nburn
           
       %Store results after burn-in    
       A_post(:,count-nburn)=A_old;
       
       AA=[1 0 0; -A_old(1,1) 1 0; -A_old(2:3,1)' 1];
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 7: Generate a draw for d(ii)^-1 from independent gamma %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       omega=AA*S*AA';
       d11=inv(gamrnd(kappastar,1/gettau(kappa,omega(1,1),AA(1,:),Syy,Sxy,eta,M1,M_star1)));
       d22=inv(gamrnd(kappastar,1/gettau(kappa,omega(2,2),AA(2,:),Syy,Sxy,eta,M2,M_star2)));
       d33=inv(gamrnd(kappastar,1/gettau(kappa,omega(3,3),AA(3,:),Syy,Sxy,eta,M3,M_star3)));
       DD=diag([d11;d22;d33]);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 8: Generate a draw for b(i) from multivariate normal %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       m1_star=getmstar(M_star1,Sxy,M1,eta,AA(1,:));
       m2_star=getmstar(M_star2,Sxy,M2,eta,AA(2,:));
       m3_star=getmstar(M_star3,Sxy,M3,eta,AA(3,:));
       
       b1=m1_star+(randn(1,nlags*n+1)*chol(d11.*M_star1))';
       b2=m2_star+(randn(1,nlags*n+1)*chol(d22.*M_star2))';
       b3=m3_star+(randn(1,nlags*n+1)*chol(d33.*M_star3))';
       BB=[b1';b2';b3']; 
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 9: Compute IRFs (not cumulated) %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       IRF(:,:,:,count-nburn)=impulse_response_1SD(AA,DD,BB,n,nlags,hmax-1);
       
       clear AA DD BB
       
    end 
    
end

% Prepare histograms 
nbin=500;

for jc=1:size(A_post,1)
    [ag,bg]=hist(A_post(jc,:),nbin);
    delta=bg(1,2)-bg(1,1);
    bg_i(jc,:)=bg;
    post_i(jc,:)=ag./((ndraws-nburn)*delta);
    clear ag bg delta
end

%compute posterior density for demand elasticity
jd=2;
[agi,bgi]=hist((1./A_post(jd,:)),5000000);
deltai=bgi(1,2)-bgi(1,1);
post_ii=agi./((ndraws-nburn)*deltai);

% Compute acceptance ratio of RW-MH algorithm
acceptance_ratio=naccept/ndraws;
disp(['Acceptance ratio:' num2str(acceptance_ratio)])

figure(2)
subplot(2,2,1)
bar(bg_i(1,:),post_i(1,:)), hold on, plot(z1,prior_eta1,'r','linewidth',2); box on
axis([-1 1 0 3.5])
title('${\alpha}_{yp}$','interpreter','latex','fontsize',18)

subplot(2,2,2)
bar(bg_i(2,:),post_i(2,:)), hold on, plot(z1,prior_beta,'r','linewidth',2); box on
axis([-1 1 0 2.5])
title('${\alpha}_{pq}$','interpreter','latex','fontsize',18)

subplot(2,2,3)
bar(bg_i(3,:),post_i(3,:)), hold on, plot(z1,prior_gamma2,'r','linewidth',2); box on
axis([0 0.5 0 7])
title('${\alpha}_{py}$','interpreter','latex','fontsize',18)

subplot(2,2,4)
bar(bgi,post_ii)
axis([-30 30 0 0.25])
title('Short-run oil demand elasticity','interpreter','latex','fontsize',13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM=mean(A_post,2);
AA=[1 0 0; -AM(1,1) 1 0; -AM(2:3,1)' 1];
HH=inv(AA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot impulse responses for one-standard deviation shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=0.025;
index=[alpha*(ndraws-nburn) (1-alpha)*(ndraws-nburn)];  %implies 68% coverage of the entire distribution
HO=(0:1:hmax-1)';    %impulse response horizon

%plots IRFs together with IRFs from Kilian AER 
%obtained with standard Cholesky decomposition
figure_IRFs_together

% figure(1)
% subplot(3,3,1)
% x1=-sort(cumsum(squeeze(IRF(1,1,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -2.5 1.5])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-2.5:1:1.5)
% ylabel('Oil production','fontsize',14)
% title('Oil supply shock','fontsize',14)
% %ylabel('percent','fontsize',12)
% 
% subplot(3,3,4)
% x1=sort(cumsum(squeeze(IRF(1,2,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -2.5 1.5])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-2.5:1:1.5)
% ylabel('Oil production','fontsize',14)
% title('Aggregate demand shock','fontsize',14)
% 
% subplot(3,3,7)
% x1=sort(cumsum(squeeze(IRF(1,3,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -2.5 1.5])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-2.5:1:1.5)
% ylabel('Oil production','fontsize',14)
% title('Oil-specific demand shock','fontsize',14)
% xlabel('Months','fontsize',12)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,2)
% x1=-sort(squeeze(IRF(2,1,:,:)),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -5 10])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-5:5:10)
% ylabel('Real activity','fontsize',14)
% title('Oil supply shock','fontsize',14)
% %ylabel('percent','fontsize',12)
% 
% subplot(3,3,5)
% x1=sort(squeeze(IRF(2,2,:,:)),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:') 
% axis([0 hmax-1 -5 10])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-5:5:10)
% ylabel('Real activity','fontsize',14)
% title('Aggregate demand shock','fontsize',14)
% 
% subplot(3,3,8)
% x1=sort(squeeze(IRF(2,3,:,:)),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:') 
% axis([0 hmax-1 -5 10])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-5:5:35)
% ylabel('Real activity','fontsize',14)
% title('Oil-specific demand shock','fontsize',14)
% xlabel('Months','fontsize',12)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(3,3,3)
% x1=-sort(squeeze(IRF(3,1,:,:)),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -5 10])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-5:5:10)
% ylabel('Real oil price','fontsize',14)
% title('Oil supply shock','fontsize',14)
% %ylabel('percent','fontsize',12)
% 
% subplot(3,3,6)
% x1=sort(squeeze(IRF(3,2,:,:)),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -5 10])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-5:5:10)
% ylabel('Real oil price','fontsize',14)
% title('Aggregate demand shock','fontsize',14)
% 
% subplot(3,3,9)
% x1=sort(squeeze(IRF(3,3,:,:)),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -5 10])
% set(gca,'XTick',0:5:20)
% set(gca,'YTick',-5:5:10)
% title('Oil-specific demand shock','fontsize',14)
% ylabel('Real oil price','fontsize',14)
% xlabel('Months','fontsize',12)
% 
% save results_KAER
