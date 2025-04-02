%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             QUESTION 2,5            %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all; close all; clc
load('data2.mat');

% Variables
gdp_growth = gdp;
inflation  = cpi;
fedfunds   = dff;
M2_growth  = m2;

% Time variables
tstart       = 1;                % start estimation with 1959:Q2
date_start   = 1985;
date_end     = 2006.12;
time         = (1984.75:0.25:2007.1)';
correct_size = size(time);
time         = time(1+1:correct_size-1);
size(time);

% Defining matrix of variables
data=[gdp_growth inflation fedfunds M2_growth ];
                   
const = ones(size(data,1),1);   % deterministic term: constant
ndet  = 1;                       % number of deterministic variables
data  = [const data];

% Model setup
[sp,nvars] = size(data);        %sp is the total number of observations
nvars      = nvars-ndet;        %the '-ndet' takes out the counting of deterministic variables
nlags      = 4;                 %number of lags

ess  = sp-nlags;        %effective sample size after taking lags into account
sb   = nlags+1;         %sample beginning
sl   = sp;              %last period of sample
ncoe = nvars*nlags;     %number of coefficients without deterministic variables a los que formular priors

% Building X for Y=X*B+U
x=data(:,ndet+1:end);
X=zeros(ess,ncoe);
for k=1:nlags
    X(:,nvars*(k-1)+1:nvars*k) = x(sb-k:sl-k,:);
end
X=[ones(ess,1) X];   %without trend

% OLS estimation
y    = x(sb:sl,:);
xtx  = X'*X;
xty  = X'*y;
% OLS estimation of reduced-form coefficients
beta = inv(X'*X)*X'*y;
e    = y-X*beta;

% Creating an empty string cell array of size 49x1
list    = cell(nlags*4,1);
list{1} = 'Constant';
k = 2;
for i = 1:nlags
    list{k} = sprintf('GDP lag %d', i);
    k = k+1;
    list{k} = sprintf('Inflation lag %d', i);
    k = k+1;
    list{k} = sprintf('Fed funds rate %d', i);
    k = k+1;
    list{k} = sprintf('M2 lag %d', i);
    k = k+1;
end
betas = array2table(beta, 'RowNames', list , 'VariableNames', {'Deflactor', 'GDP', 'Fed funds', 'M2'});

% Variance-covariance matrix 
vmat  = e'*e/ess;
covariance_matrix = array2table(vmat, 'RowNames', {'Deflactor', 'GDP', 'Fed funds', 'M2'} , 'VariableNames', {'Deflactor', 'GDP', 'Fed funds', 'M2'});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%            ALGORITHM FOR GETTING POSTERIORS           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ndraws=2e5;     %number of Metropolis-Hastings iterations
nburn=1e5;      %number of burn-in draws

xsi=0.75^2;       %tuning parameter
hmax=20;          %impulse response horizon (20 quarters // 5 years)

yy=[gdp_growth inflation fedfunds M2_growth ]; %Rearrange data
s=size(time,1); % Observations per variables
n=size(yy,2); %Number of variables

%Data matrices
[X,y,T]=getXy(yy,ndet,nlags);
yyy=y;
xxx=X;

%Set seed
seednumber=07081924;
rand('seed',seednumber);
randn('seed',seednumber);

%%% BULDING PRIORS %%% 

% For matrix A
bounds = 5;
z1=-bounds:.01:bounds; %very no restrictive grid for parameters where no sign is imposed a priori
% We impose the same prior for all the variables below
c_general = 0;
sigma_general = 100;   
nu_general = 3;
% Setting priors:
prior_beta1  = student_prior(z1,c_general,sigma_general,nu_general);
prior_theta1 = student_prior(z1,c_general,sigma_general,nu_general);
prior_theta2 = student_prior(z1,c_general,sigma_general,nu_general);
prior_gamma1 = student_prior(z1,c_general,sigma_general,nu_general);
prior_gamma2 = student_prior(z1,c_general,sigma_general,nu_general);
prior_gamma3 = student_prior(z1,c_general,sigma_general,nu_general);
% Set arbitrary initial values for elements of A 
A_old=[c_general; c_general; c_general; c_general;c_general;c_general];   %prior mode of elements in A
c=size(A_old,1);                    %number of parameters to be estimated in A
% The full matrix A for construct is
% A =[1 0 0 0; -A_old(1,1) 1 0 0; -A_old(2:3,1)' 1 0; -A_old(4:6,1)' 1 ];

% For lagged coefficients B

% Compute standard deviation of each series residual via an OLS regression
% to be used in setting the prior (here: AR(4))
[s11,uhat1] = sd_prior(yy(:,1),nlags);
[s22,uhat2] = sd_prior(yy(:,2),nlags);
[s33,uhat3] = sd_prior(yy(:,3),nlags);
[s44,uhat4] = sd_prior(yy(:,4),nlags);

%Selecting hyperparameters for our Minnesota prior (and to make it uninforative)
lambda0 = 10^9;    %overall confidence in prior 
lambda1 = 1;       %confidence on higher-order lags
lambda2 = 1;       %confidence in other-than-own lags 
lambda3 = 100;     %tightness of constant term 
% Prior mean (a RW Minnesota prior)
eta = [eye(n) zeros(n,(n*nlags+1)-n)];
% Prior covariance
SS  = [s11;s22;s33;s44];
M   = getM(lambda0,lambda1,lambda2,lambda3,SS,nlags,n);
M1  = M(:,:,1);
M2  = M(:,:,2);
M3  = M(:,:,3);
M4  = M(:,:,4);

% Compute summary statistics of the observed data:
Syy=yyy'*yyy;
Sxx=xxx'*xxx;
Sxy=yyy'*xxx;
omega_hatT=(Syy-Sxy*inv(Sxx)*Sxy')/T;
% Compute M_star(i) (posterior covariance) 
M_star1=inv(Sxx+inv(M1));
M_star2=inv(Sxx+inv(M2));
M_star3=inv(Sxx+inv(M3));
M_star4=inv(Sxx+inv(M4));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1c: Set priors for inverse of diagonal elements (D) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% The mean is calibrated on diagonal elements in omega
uhat=[uhat1 uhat2 uhat3 uhat4];
S=uhat'*uhat/T;
kappa=0.5; % a very small parameter for an uninformative prior
kappastar=kappa+(T/2); % posterior kappastar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hamilton's proposal to deal with scaling issue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixed parameter values
param  = [c_general;sigma_general;nu_general;c_general;sigma_general;nu_general; c_general;sigma_general;nu_general; ...
c_general;sigma_general;nu_general;c_general;sigma_general;nu_general;c_general;sigma_general;nu_general;T;vec(omega_hatT);kappastar;kappa];   
f_anon = @(theta_hat)post_val(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M4,M_star1,M_star2,M_star3,M_star4,S);

% find posterior mode
theta_zero = A_old;     %mode of prior distributions
options = optimset('LargeScale','off'); %optimization options
[theta_max,val_max,exitm,om,gm,HM] = fminunc(f_anon,theta_zero,options);
PH=chol(inv(HM))';

%start MH algorithm with theta_hat
A_old=theta_max;
A=[1 0 0 0; -A_old(1,1) 1 0 0; -A_old(2:3,1)' 1 0; -A_old(4:6,1)' 1 ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Set the variance of the candidate generating density (P) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the variance of the candidate generating density (P)
W=xsi*eye(c);     %variance of RW-MH  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Evaluate posterior at starting value for A:  % (pag 104)(we compute all inputs necesary to evaluate the posterior in A)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

omega=A*S*A';
taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);
% Evaluate prior p(A) at old draw
prior_beta1_eval = student_prior(A_old(1,1),c_general,sigma_general,nu_general);
prior_theta1_eval = student_prior(A_old(2,1),c_general,sigma_general,nu_general);
prior_theta2_eval = student_prior(A_old(3,1),c_general,sigma_general,nu_general);
prior_gamma1_eval = student_prior(A_old(4,1),c_general,sigma_general,nu_general);
prior_gamma2_eval = student_prior(A_old(5,1),c_general,sigma_general,nu_general);
prior_gamma3_eval = student_prior(A_old(6,1),c_general,sigma_general,nu_general);
% Compute posterior value at candidate draw
log_priors=log(prior_beta1_eval)+log(prior_theta1_eval)+log(prior_theta2_eval)+log(prior_gamma1_eval)+log(prior_gamma2_eval)+log(prior_gamma3_eval);
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
       prior_beta1_eval = student_prior(A_new(1,1),c_general,sigma_general,nu_general);
       prior_theta1_eval = student_prior(A_new(2,1),c_general,sigma_general,nu_general);
       prior_theta2_eval = student_prior(A_new(3,1),c_general,sigma_general,nu_general);
       prior_gamma1_eval = student_prior(A_new(4,1),c_general,sigma_general,nu_general);
       prior_gamma2_eval = student_prior(A_new(5,1),c_general,sigma_general,nu_general);
       prior_gamma3_eval = student_prior(A_new(6,1),c_general,sigma_general,nu_general);
       
       % Construct full matrix A 
       A=[1 0 0 0; -A_new(1,1) 1 0 0; -A_new(2:3,1)' 1 0; -A_new(4:6,1)' 1 ];

       omega=A*S*A';
       taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
       taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
       taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
       taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);
       % Compute posterior value at new candidate draw
       log_priors=log(prior_beta1_eval)+log(prior_theta1_eval)+log(prior_theta2_eval)+log(prior_gamma1_eval)+log(prior_gamma2_eval)+log(prior_gamma3_eval);
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
       
       AA=[1 0 0 0; -A_old(1,1) 1 0 0; -A_old(2:3,1)' 1 0; -A_old(4:6,1)' 1];
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 7: Generate a draw for d(ii)^-1 from independent gamma %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       omega=AA*S*AA';
       d11=inv(gamrnd(kappastar,1/gettau(kappa,omega(1,1),AA(1,:),Syy,Sxy,eta,M1,M_star1)));
       d22=inv(gamrnd(kappastar,1/gettau(kappa,omega(2,2),AA(2,:),Syy,Sxy,eta,M2,M_star2)));
       d33=inv(gamrnd(kappastar,1/gettau(kappa,omega(3,3),AA(3,:),Syy,Sxy,eta,M3,M_star3)));
       d44=inv(gamrnd(kappastar,1/gettau(kappa,omega(4,4),AA(4,:),Syy,Sxy,eta,M4,M_star4)));       
       DD=diag([d11;d22;d33;d44]);
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % STEP 8: Generate a draw for b(i) from multivariate normal %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       m1_star=getmstar(M_star1,Sxy,M1,eta,AA(1,:));
       m2_star=getmstar(M_star2,Sxy,M2,eta,AA(2,:));
       m3_star=getmstar(M_star3,Sxy,M3,eta,AA(3,:));
       m4_star=getmstar(M_star4,Sxy,M4,eta,AA(4,:));
       
       b1=m1_star+(randn(1,nlags*n+1)*chol(d11.*M_star1))';
       b2=m2_star+(randn(1,nlags*n+1)*chol(d22.*M_star2))';
       b3=m3_star+(randn(1,nlags*n+1)*chol(d33.*M_star3))';
       b4=m4_star+(randn(1,nlags*n+1)*chol(d44.*M_star4))';
       BB=[b1';b2';b3';b4']; 
       
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
subplot(2,3,1)
bar(bg_i(1,:),post_i(1,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_beta1,'r','linewidth',2); box on
line([median(bg_i(1,:)) median(bg_i(1,:))], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
xlim([-0.5 0.2])
title('${\beta}_{1}$','interpreter','latex','fontsize',28)

subplot(2,3,2)
bar(bg_i(2,:),post_i(2,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_theta1,'r','linewidth',2); box on
line([median(bg_i(2,:)) median(bg_i(2,:))], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
xlim([-0.25 0.25])
title('${\theta}_{1}$','interpreter','latex','fontsize',28)

subplot(2,3,3)
bar(bg_i(3,:),post_i(3,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_theta2,'r','linewidth',2); box on
line([median(bg_i(3,:)) median(bg_i(3,:))], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
xlim([0 0.6])
title('${\theta}_{2}$','interpreter','latex','fontsize',28)
legend('Posterior distribution','Prior distribution')

subplot(2,3,4)
bar(bg_i(4,:),post_i(4,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_gamma1,'r','linewidth',2); box on
line([median(bg_i(4,:)) median(bg_i(4,:))], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
xlim([-0.3 0.45])
title('${\gamma}_{1}$','interpreter','latex','fontsize',28)

subplot(2,3,5)
bar(bg_i(5,:),post_i(5,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_gamma2,'r','linewidth',2); box on
line([median(bg_i(5,:)) median(bg_i(5,:))], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
xlim([-1.7 -0.7])
title('${\gamma}_{2}$','interpreter','latex','fontsize',28)

subplot(2,3,6)
bar(bg_i(6,:),post_i(6,:),'FaceColor', [0.000000 0.200000 0.600000]), hold on, plot(z1,prior_gamma3,'r','linewidth',2); box on
line([median(bg_i(6,:)) median(bg_i(6,:))], ylim, 'Color', 'g', 'LineStyle', '--', 'LineWidth', 2);
xlim([-1.3 0.1])
title('${\gamma}_{3}$','interpreter','latex','fontsize',28)

set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1400, 900]);

if ~exist('graphs 3', 'dir')
    mkdir('graphs 3');
end

saveas(gcf, 'graphs 3/Posterior distribution question 5 .png');


%Because all the coefficients have the same prior, we only graph one of them
figure(3)
plot(z1,prior_gamma2,'r','linewidth',2);
title('Priors for all the coefficients','interpreter','latex','fontsize',13)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AM=mean(A_post,2);
AA=[1 0 0 0; -AM(1,1) 1 0 0; -AM(2:3,1)' 1 0; -AM(4:6,1)' 1];
HH=inv(AA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot impulse responses for one-standard deviation shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alpha=0.16;
index=[alpha*(ndraws-nburn) (1-alpha)*(ndraws-nburn)];  %implies 68% coverage of the entire distribution
HO=(0:1:hmax-1)';    %impulse response horizon

%%%%%%%%%%%%%%%%%%%%%%%% Response of US GDP %%%%%%%%%%%%%%%%%%%%%%%%%%
figure(4)
subplot(4,4,1)
x1=-sort(cumsum(squeeze(-IRF(1,1,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('GDP','fontsize',14)
title('GDP shock','fontsize',14)

subplot(4,4,5)
x1=sort(cumsum(squeeze(-IRF(1,2,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('GDP','fontsize',14)
title('Inflation shock','fontsize',14)

subplot(4,4,9)
x1=sort(cumsum(squeeze(-IRF(1,3,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('GDP','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,13)
x1=sort(cumsum(squeeze(IRF(1,4,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('GDP','fontsize',14)
title('M2 shock','fontsize',14)


%%%%%%%%%%%%%%%%%%%%%%%% Inflation Response %%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,4,2)
x1=-sort(cumsum(squeeze(-IRF(2,1,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Inflation','fontsize',14)
title('GDP shock','fontsize',14)

subplot(4,4,6)
x1=sort(cumsum(squeeze(-IRF(2,2,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Inflation','fontsize',14)
title('Inflation shock','fontsize',14)

subplot(4,4,10)
x1=sort(cumsum(squeeze(-IRF(2,3,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Inflation','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,14)
x1=sort(cumsum(squeeze(IRF(2,4,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis([0 hmax-1 -2.5 1.5])
ylabel('Inflation','fontsize',14)
title('M2 shock','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Response of interest rate %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,4,3)
x1=-sort(squeeze(-IRF(3,1,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Fed rate','fontsize',14)
title('GDP shock','fontsize',14)

subplot(4,4,7)
x1=-sort(squeeze(-IRF(3,2,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Fed rate','fontsize',14)
title('Inflation shock','fontsize',14)

subplot(4,4,11)
x1=-sort(squeeze(-IRF(3,3,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Fed rate','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,15)
x1=-sort(squeeze(IRF(3,4,:,:)),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('Fed rate','fontsize',14)
title('M2 shock','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Response of Money Supply %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,4,4)
x1=-sort(cumsum(squeeze(-IRF(4,1,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('M2','fontsize',14)
title('GDP shock','fontsize',14)

subplot(4,4,8)
x1=sort(cumsum(squeeze(-IRF(4,2,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('M2','fontsize',14)
title('Inflation shock','fontsize',14)

subplot(4,4,12)
x1=sort(cumsum(squeeze(-IRF(4,3,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('M2','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,16)
x1=sort(cumsum(squeeze(IRF(4,4,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
axis tight
ylabel('M2','fontsize',14)
title('M2 shock','fontsize',14)
set(gcf, 'PaperPositionMode', 'auto');
set(gcf, 'Position', [100, 100, 1400, 900]);
if ~exist('graphs 3', 'dir')
    mkdir('graphs 3');
end

% Saving as .png
saveas(gcf, 'graphs 3/question 5 IRF .png');