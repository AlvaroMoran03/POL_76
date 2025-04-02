
clear all;
clc;

% Downloading the data from the FRED Web Page
V1 = getFredData('GDPC1', '1960-01-01', '2006-12-31','pch','q');
V2 = getFredData('GDPDEF', '1960-01-01', '2006-12-31','pch','q');
V3 = getFredData('M2REAL', '1960-01-01', '2006-12-31','pch','q');
V4 = getFredData('FEDFUNDS', '1960-01-01', '2006-12-31','lin','q','avg');

% Renaming the Variables
GDP = V1.Data(:,2); DEF = V2.Data(:,2); M2 = V3.Data(:,2); FED = V4.Data(:,2);

clear V1 V2 V3 V4

% Setting a matrix that contains requiered time series
Data=[GDP DEF FED M2];

time = (1960:0.25:2006.75)';

yy = Data;
clearvars -except yy time

ndraws=1e6;     % number of MH iterations 2e5
nburn=5e5;      % number of burn-in draws 1e5

% NOTE:
% The parameter xsi is the tuning parameter for the MH step and has to be
% re-set for each application such that it implies a 30-35% acceptance ratio.
% If acceptance ratio is too high, make xsi bigger; 
% if acceptance ratio is too low, make xsi smaller.

xsi=0.7^2;       % tuning parameter to achieve a target acceptance rate of around 35%

nlags=4;         % number of lags
hmax=20;         % impulse response horizon (17 months)
ndet=1;          % number of deterministic variables 
                 % (1: constant; 2: constant and time trend)

s=size(time,1);
n=size(yy,2);                 % number of endogenous variables
                              
% Get data matrices
[X,y,T]=getXy(yy,ndet,nlags); % input: yy (las series de tiempo)
yyy=y; 
xxx=X; % la matriz de rezagos

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
z1=-bounds:.01:bounds;   % grid for parameters where no sign is imposed a priori

% Setting priors:

alpha1 = 0;
sigma_alpha1 = 100;   
nu_alpha1 = 3;
prior_alpha1 = student_prior(z1,alpha1,sigma_alpha1,nu_alpha1);

beta1 = 0;
sigma_beta1 = 100;   
nu_beta1 = 3;
prior_beta1 = student_prior(z1,beta1,sigma_beta1,nu_beta1);

beta2 = 0;
sigma_beta2 = 100;   
nu_beta2 = 3;
prior_beta2 = student_prior(z1,beta2,sigma_beta2,nu_beta2);

delta1 = 0;
sigma_delta1 = 100;   
nu_delta1 = 3;
prior_delta1 = student_prior(z1,delta1,sigma_delta1,nu_delta1);

delta2 = 0;
sigma_delta2 = 100;   
nu_delta2 = 3;
prior_delta2 = student_prior(z1,delta2,sigma_delta2,nu_delta2);

delta3 = 0;
sigma_delta3 = 100;   
nu_delta3 = 3;
prior_delta3 = student_prior(z1,delta3,sigma_delta3,nu_delta3);

% Set arbitrary initial values for elements of A 

A_old = [alpha1; beta1; beta2; delta1; delta2; delta3]; % prior mode of elements in A
c=size(A_old,1);  % number of parameters to be estimated in A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1b: Set informative priors on lagged coefficients (B) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute standard deviation of each series residual via an OLS regression

[s11,uhat1]=sd_prior(yy(:,1),nlags); % sii: desv estandar muestral del error de la variable ii
[s22,uhat2]=sd_prior(yy(:,2),nlags);
[s33,uhat3]=sd_prior(yy(:,3),nlags);
[s44,uhat4]=sd_prior(yy(:,4),nlags);

% See Doan (2013) for choices of values of hyperparameters (Minnesota prior)
lambda0=10^9;    % overall confidence in prior (smaller lambda0 implies greater weight to RW) (little confid in prior bec the orig study dont use a bay approach, only rely on data)
lambda1=1;       % confidence on higher-order lags (lambda1 = 0 gives all lags equal weight)
lambda2=1;       % confidence in other-than-own lags 
lambda3=100;     % tightness of constant term 


% EXPECTED VALUE OF REDUCED-FORM COEFFICIENTS (page 91): RW prior
eta=[eye(n) zeros(n,(n*nlags+1)-n)]; % doesn't matter if we set the growth rate of production as a RW bec we give poca imp a prior)

% PRIOR COVARIANCE
SS=[s11;s22;s33;s44]; 
M=getM(lambda0,lambda1,lambda2,lambda3,SS,nlags,n); 
M1=M(:,:,1); % var component of "i" equation
M2=M(:,:,2);
M3=M(:,:,3);
M4=M(:,:,4);

% Compute summary statistics of the observed data:
Syy=yyy'*yyy;
Sxx=xxx'*xxx;
Sxy=yyy'*xxx;
omega_hatT=(Syy-Sxy*inv(Sxx)*Sxy')/T;

% Compute M_star(i) (posterior variance of B) (pag 94)
M_star1=inv(Sxx+inv(M1)); % we combine info in data with prior to get M_star 
M_star2=inv(Sxx+inv(M2));
M_star3=inv(Sxx+inv(M3));
M_star4=inv(Sxx+inv(M4));
% M_star1=inv(Sxx);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1c: Set priors for inverse of diagonal elements (D) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The mean is calibrated on diagonal elements in omega
uhat=[uhat1 uhat2 uhat3 uhat4];
S=uhat'*uhat/T; % matriz de var-cov de los residuos del VAR reducido

kappa=0.5; % bajo kappa --> le estamos dando menos importancia a prior, es como tener, solo 1 obs
kappastar=kappa+(T/2); % posterior kappastar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JIM's PROPOSAL to deal with scaling issue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fixed parameter values

param= [alpha1;sigma_alpha1;nu_alpha1;beta1;sigma_beta1;nu_beta1; ...
         beta2;sigma_beta2;nu_beta2;delta1;sigma_delta1;nu_delta1; ... 
         delta2;sigma_delta2;nu_delta2;delta3;sigma_delta3;nu_delta3; ...
         T;vec(omega_hatT);kappastar;kappa];

f_anon = @(theta_hat)post_val(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M4,M_star1,M_star2,M_star3,M_star4,S);

%find posterior mode

theta_zero = A_old;     % mode of prior distributions
options = optimset('LargeScale','off');
[theta_max,val_max,exitm,om,gm,HM] = fminunc(f_anon,theta_zero,options);
PH=chol(inv(HM))';

% start MH algorithm with theta_hat
A_old=theta_max;
%A = setA(A_old);
A = [1 0 0 0; -A_old(1,1) 1 0 0; -A_old(2:3,1)' 1 0; -A_old(4:6,1)' 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Set the variance of the candidate generating density (P) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

W=xsi*eye(c);     % variance of RW-MH  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Evaluate posterior at starting value for A:  % (pag 104)(we compute all inputs necesary to evaluate the posterior in A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=A*S*A';
taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);

% Evaluate prior p(A) at old draw

prior_a = student_prior(A_old(1,1),alpha1,sigma_alpha1,nu_alpha1);
prior_b = student_prior(A_old(2,1),beta1,sigma_beta1,nu_beta1);
prior_c = student_prior(A_old(3,1),beta2,sigma_beta2,nu_beta2);
prior_d = student_prior(A_old(4,1),delta1,sigma_delta1,nu_delta1);
prior_e = student_prior(A_old(5,1),delta2,sigma_delta2,nu_delta2);
prior_f = student_prior(A_old(6,1),delta3,sigma_delta3,nu_delta3);

% Compute posterior value at candidate draw
log_priors = log(prior_a) + log(prior_b) + log(prior_c) ...
            + log(prior_d) + log(prior_e) + log(prior_f);
up=log_priors+T/2*log(det(A*omega_hatT*A'));
down=kappastar*log((2/T)*taustar1) ...
    +kappastar*log((2/T)*taustar2) ...
    +kappastar*log((2/T)*taustar3) ...
    +kappastar*log((2/T)*taustar4);
posteriorOLD=up-down;

% RW-MH algorithm 
naccept=0;
count=0;

% Store posterior distribution (after burn-in)
A_post=zeros(c,ndraws-nburn);                    % remember (c = # de params a ser estimados en A)
IRF=zeros(n,n,hmax,ndraws-nburn);

while count<ndraws 
      count=count+1;
      if (count/10000) == floor(count/10000)
          count
      end
      
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
    % STEP 4a: Generate draw for A from the RW candidate density %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    A_new=A_old+chol(W)'*PH*randn(c,1)/sqrt(0.5*(randn(1)^2 + randn(1)^2));  % fat tails
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % STEP 4b: Evaluate posterior at new draw %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

       % Evaluate prior p(A) at new draw
       prior_a = student_prior(A_new(1,1),alpha1,sigma_alpha1,nu_alpha1);
       prior_b = student_prior(A_new(2,1),beta1,sigma_beta1,nu_beta1);
       prior_c = student_prior(A_new(3,1),beta2,sigma_beta2,nu_beta2);
       prior_d = student_prior(A_new(4,1),delta1,sigma_delta1,nu_delta1);
       prior_e = student_prior(A_new(5,1),delta2,sigma_delta2,nu_delta2);
       prior_f = student_prior(A_new(6,1),delta3,sigma_delta3,nu_delta3);

       % Construct full matrix A 
       %A = setA(A_new);
       A = [1 0 0 0; -A_new(1,1) 1 0 0; -A_new(2:3,1)' 1 0; -A_new(4:6,1)' 1];

       omega=A*S*A';
       taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
       taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
       taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
       taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);
           
       % Compute posterior value at new candidate draw
       log_priors = log(prior_a) + log(prior_b) + log(prior_c) ...
                    + log(prior_d) + log(prior_e) + log(prior_f);
       
       up=log_priors+T/2*log(det(A*omega_hatT*A'));
       down=kappastar*log((2/T)*taustar1) ...
           +kappastar*log((2/T)*taustar2) ...
           +kappastar*log((2/T)*taustar3) ...
           +kappastar*log((2/T)*taustar4);
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
       %AA= setA(A_post);
       AA = [1 0 0 0; -A_old(1,1) 1 0 0; -A_old(2:3,1)' 1 0; -A_old(4:6,1)' 1];

       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 7: Generate a draw for d(ii)^-1 from independent gamma % (page 95)
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       omega=AA*S*AA';
       d11=inv(gamrnd(kappastar,1/gettau(kappa,omega(1,1),AA(1,:),Syy,Sxy,eta,M1,M_star1))); % matlab function computes 1/tau, therefore: 1/gettau is tau*
       d22=inv(gamrnd(kappastar,1/gettau(kappa,omega(2,2),AA(2,:),Syy,Sxy,eta,M2,M_star2)));
       d33=inv(gamrnd(kappastar,1/gettau(kappa,omega(3,3),AA(3,:),Syy,Sxy,eta,M3,M_star3)));
       d44=inv(gamrnd(kappastar,1/gettau(kappa,omega(4,3),AA(4,:),Syy,Sxy,eta,M4,M_star4)));
       DD=diag([d11;d22;d33;d44]);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 8: Generate a draw for b(i) from multivariate normal % (page 94)
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       

       m1_star=getmstar(M_star1,Sxy,M1,eta,AA(1,:)); %cambiar por estima OLS
       m2_star=getmstar(M_star2,Sxy,M2,eta,AA(2,:));
       m3_star=getmstar(M_star3,Sxy,M3,eta,AA(3,:));
       m4_star=getmstar(M_star4,Sxy,M4,eta,AA(4,:));

       % b_ols = inv(xxx'*xxx)*(xxx'*yyy);
       % m1_star=b_ols(:,1);
       % m2_star=b_ols(:,2);
       % m3_star=b_ols(:,3);
       % m4_star=b_ols(:,4);

       b1=m1_star+(randn(1,nlags*n+1)*chol(d11.*M_star1))'; % bi = mi_star + N (0,diiM*i) -> pq le aplica cholesky?
       b2=m2_star+(randn(1,nlags*n+1)*chol(d22.*M_star2))';
       b3=m3_star+(randn(1,nlags*n+1)*chol(d33.*M_star3))';
       b4=m4_star+(randn(1,nlags*n+1)*chol(d44.*M_star4))';
       BB=[b1';b2';b3';b4']; 
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 9: Compute IRFs (not cumulated) %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
       IRF(:,:,:,count-nburn)=impulse_response_1SD(AA,DD,BB,n,nlags,hmax-1); % for a draw of A, D, B we take a transformation
       
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
[agi,bgi]=hist((1./A_post(jd,:)),5000000); % comput demand elasticity as the inverse of each element
deltai=bgi(1,2)-bgi(1,1);
post_ii=agi./((ndraws-nburn)*deltai);

% Compute acceptance ratio of RW-MH algorithm
acceptance_ratio=naccept/ndraws;
disp(['Acceptance ratio:' num2str(acceptance_ratio)])


figure(2)
subplot(2,3,1)
bar(bg_i(1,:),post_i(1,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_alpha1,'r','linewidth',2); box on
xlim([-0.2 0.2])
title('${\alpha}_{1}$','interpreter','latex','fontsize',28)

subplot(2,3,2)
bar(bg_i(2,:),post_i(2,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_beta1,'r','linewidth',2); box on
xlim([-0.2 0.5])
title('${\beta}_{1}$','interpreter','latex','fontsize',28)

subplot(2,3,3)
bar(bg_i(3,:),post_i(3,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_beta2,'r','linewidth',2); box on
xlim([-0.2 1.8])
title('${\beta}_{2}$','interpreter','latex','fontsize',28)
legend('Posterior distribution','Prior distribution')

subplot(2,3,4)
bar(bg_i(4,:),post_i(4,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_delta1,'r','linewidth',2); box on
xlim([-0.2 0.35])
title('${\delta}_{1}$','interpreter','latex','fontsize',28)

subplot(2,3,5)
bar(bg_i(5,:),post_i(5,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_delta2,'r','linewidth',2); box on
xlim([-1.5 0.2])
title('${\delta}_{2}$','interpreter','latex','fontsize',28)

subplot(2,3,6)
bar(bg_i(6,:),post_i(6,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_delta3,'r','linewidth',2); box on
xlim([-0.6 -0])
title('${\delta}_{3}$','interpreter','latex','fontsize',28)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

AM=mean(A_post,2);
AA=setA(AM);
HH=inv(AA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot impulse responses for one-standard deviation shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

alpha=0.025;
index=[alpha*(ndraws-nburn) (1-alpha)*(ndraws-nburn)];  % implies 68% coverage of the entire distribution
HO=(0:1:hmax-1)';    %impulse response horizon

% Plot IRFs, obtained with standard Cholesky decomposition

figure(3)

%%%%%%%%%%%%%%%%%%%%%%%% Response of US GDP %%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,4,1)
x1= sort(squeeze(IRF(1,1,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.8])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:0.5:1.5)
ylabel('GDP','fontsize',16)
title('Demand shock','fontsize',38)
ax = gca;
ax.FontSize = 14;

subplot(4,4,5)
x1= sort(squeeze(IRF(1,2,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.8])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:0.5:1.5)
ylabel('GDP','fontsize',16)
title('Supply shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,9)
x1= sort(squeeze(IRF(1,3,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.8])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:0.5:1.5)
ylabel('GDP','fontsize',16)
%title('M2 shock','fontsize',11)
title('MP shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,13)
x1= sort(squeeze(IRF(1,4,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.8])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.5:0.5:1.5)
ylabel('GDP','fontsize',16)
%title('MP Shock','fontsize',11)
title('Liquidity shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%% Inflation Response %%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,4,2)
x1= sort(squeeze(IRF(2,1,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.2 0.3])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.2:2.0)
ylabel('Inflation','fontsize',16)
title('Demand Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,6)
x1= sort(squeeze(IRF(2,2,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.2 0.3])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.2:2.0)
ylabel('Inflation','fontsize',16)
title('Supply Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,10)
x1= sort(squeeze(IRF(2,3,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.2 0.3])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.2:2.0)
ylabel('Inflation','fontsize',16)
%title('M2 shock','fontsize',11)
title('MP Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,14)
x1= sort(squeeze(IRF(2,4,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.2 0.3])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.2:2.0)
ylabel('Inflation','fontsize',16)
%title('MP shock','fontsize',11)
title('Liquidity Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Response of interest rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,4,3)
x1= sort(squeeze(IRF(3,1,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 1.2])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.5:2.0)
%ylabel('M2','fontsize',11)
ylabel('Interest rate','fontsize',16)
title('Demand Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,7)
x1= sort(squeeze(IRF(3,2,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 1.2])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.5:2.0)
%ylabel('M2','fontsize',11)
ylabel('Interest rate','fontsize',16)
title('Supply shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,11)
x1= sort(squeeze(IRF(3,3,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 1.2])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.5:2.0)
%ylabel('M2','fontsize',11)
ylabel('Interest rate','fontsize',16)
%title('M2 shock','fontsize',11)
title('MP shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,15)
x1= sort(squeeze(IRF(3,4,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 1.2])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-2.0:0.5:2.0)
%ylabel('M2','fontsize',11)
ylabel('Interest rate','fontsize',16)
%title('MP shock','fontsize',11)
title('Liquidity shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Response of Money Supply %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

subplot(4,4,4)
x1= sort(squeeze(IRF(4,1,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-0.5:0.25:0.5)
%ylabel('Fed Funds Rate','fontsize',11)
ylabel('M2','fontsize',16)
title('Demand Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,8)
x1= sort(squeeze(IRF(4,2,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-0.5:0.25:0.5)
%ylabel('Fed Funds Rate','fontsize',11)
ylabel('M2','fontsize',16)
title('Supply shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,12)
x1= sort(squeeze(IRF(4,3,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1),
axis([0 hmax-1 -0.5 0.5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-0.5:0.25:0.5)
%ylabel('Fed Funds Rate','fontsize',11)
ylabel('M2','fontsize',16)
%title('M2 shock','fontsize',11)
title('MP Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;

subplot(4,4,16)
x1= sort(squeeze(IRF(4,4,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:','linewidth',1)
axis([0 hmax-1 -0.5 0.5])
set(gca,'XTick',0:5:20)
set(gca,'YTick',-0.5:0.25:0.5)
%ylabel('Fed Funds Rate','fontsize',11)
ylabel('M2','fontsize',16)
%title('MP shock','fontsize',11)
title('Liquidity Shock','fontsize',18)
ax = gca;
ax.FontSize = 14;




