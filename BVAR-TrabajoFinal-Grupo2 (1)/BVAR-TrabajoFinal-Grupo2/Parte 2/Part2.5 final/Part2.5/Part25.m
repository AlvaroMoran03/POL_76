%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% QUESTION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

seednumber=372398;
rand('seed',seednumber);
randn('seed',seednumber);

%% 2.1 Collect quarterly data for US real GDP, the GDP deflator, the M2 money
    % stock, and the effective federal funds rate for the period 1959Q1-2006Q4
    clear all;
    clc;

    % a.Transform the data to quarterly growth rates where appropriate so that they have a useful economic
    %   interpretation. Plot the transformed data with appropriate labels.

    % Downloading the data from the FRED Web Page
    V1 = getFredData('GDPC1', '1959-04-01', '2006-12-31','pch','q');
    V2 = getFredData('GDPDEF', '1959-04-01', '2006-12-31','pch','q');
    V3 = getFredData('M2REAL', '1959-04-01', '2006-12-31','pch','q');
    V4 = getFredData('FEDFUNDS', '1959-04-01', '2006-12-31','lin','q','avg');

    % Renaming the Variables

    GDP = V1.Data(:,2); DEF = V2.Data(:,2); M2 = V3.Data(:,2); FED = V4.Data(:,2);
    
    clear V1 V2 V3 V4
    
    data=[GDP DEF FED M2];

%% 2.5 Suppose you wanted to identify the shocks underlying this model by means of sign restrictions

    % a. Provide a plot for the impact effect of a one-standard deviation shock using the analytical
      %expression for the implicit prior distribution.

      const=ones(size(data,1),1);   %deterministic term: constant
      ndet=1;                       %number of deterministic variables
      data=[const data];

      nsteps=1;                 %horizon for impulse responses                          
      ndraws=1000000;             %number of solutions that match the sign restrictions
      bmed=0.5;                 %impulse response percentile (median)
      bupp=0.84;                %upper error band (percentile)
      blow=0.16;                %lower error band (percentile)
           
      %model setup
      [sp,nvars]=size(data);   %sp is the total number of observations
      nvars=nvars-ndet;        %the '-ndet' takes out the counting of deterministic variables
      nlags=4;                 %number of lags

      ess=sp-nlags;       %effective sample size after taking lags into account
      sb=nlags+1;         %sample beginning
      sl=sp;              %last period of sample
      ncoe=nvars*nlags;   %number of coefficients without deterministic variables

      %construct X for Y=X*B+U
        
      x=data(:,ndet+1:end);
      X=zeros(ess,ncoe);
      for k=1:nlags
          X(:,nvars*(k-1)+1:nvars*k) = x(sb-k:sl-k,:);
      end
      X=[ones(ess,1) X];   %without trend

      %estimation
      y=x(sb:sl,:);
      xtx=X'*X;
      xty=X'*y;
      %OLS estimation of reduced-form coefficients
      %Bh=xr\(xr'\xty);     %xr\xr' is a faster alternative for inv(X'X)
      beta=inv(X'*X)*X'*y;
      e=y-X*beta;
      %variance-covariance matrix 
      vmat=e'*e/ess;

      %Monte Carlo simulation

      tel=0;
      count=0;

      disp('Starting Monte Carlo Simulations')
      while count<ndraws
        count=count+1
    
      %rotation matrix
      K=normrnd(0,1,nvars,nvars);   %draws an nvars x nvars matrix from a standard normal distribution
    
      [Q,R]=qr(K);                  %qr decomposition: Q is a unitary matrix
                                  %alternatively to qr decomp, you could
      for i=1:nvars                 %just calculate Q=orth(K) which delivers orthonormal matrix
         if R(i,i)<0
             Q(:,i)=-Q(:,i);
         end
      end
    
    a0=chol(vmat)'*Q;     %for point estimation
          
    %compute candidate impulse responses:
    irf(:,:,:,count)=impulse(beta',a0,ndet,nlags,nsteps);        %for point estimation                       
         
    end

    disp(['Number of full identifications: ' num2str(tel)])
    disp(['Total number of draws and rotations: ' num2str(count)])


    % Impact effect of one-standard-deviation shocks (histogram)
    figure(2)
    subplot(2,2,1)
    hist(squeeze(irf(1,1,1,:)),500)
    title('GDP')
    subplot(2,2,2)
    hist(squeeze(irf(2,1,1,:)),500)
    title('Inflation')
    subplot(2,2,3)
    hist(squeeze(irf(3,1,1,:)),500)
    title('FFR')
    subplot(2,2,4)
    hist(squeeze(irf(4,1,1,:)),500)
    title('M2')

    figure(3)
    subplot(2,2,1)
    hist(squeeze(irf(1,2,1,:)),500)
    title('GDP')
    subplot(2,2,2)
    hist(squeeze(irf(2,2,1,:)),500)
    title('Inflation')
    subplot(2,2,3)
    hist(squeeze(irf(3,2,1,:)),500)
    title('FFR')
    subplot(2,2,4)
    hist(squeeze(irf(4,2,1,:)),500)
    title('M2')

    figure(4)
    subplot(2,2,1)
    hist(squeeze(irf(1,3,1,:)),500)
    title('GDP')
    subplot(2,2,2)
    hist(squeeze(irf(2,3,1,:)),500)
    title('Inflation')
    subplot(2,2,3)
    hist(squeeze(irf(3,3,1,:)),500)
    title('FFR')
    subplot(2,2,4)
    hist(squeeze(irf(4,3,1,:)),500)
    title('M2')

    figure(5)
    subplot(2,2,1)
    hist(squeeze(irf(1,4,1,:)),500)
    title('GDP')
    subplot(2,2,2)
    hist(squeeze(irf(2,4,1,:)),500)
    title('Inflation')
    subplot(2,2,3)
    hist(squeeze(irf(3,4,1,:)),500)
    title('FFR')
    subplot(2,2,4)
    hist(squeeze(irf(4,4,1,:)),500)
    title('M2')
    
    
    vmat


    %% b.Verify empirically what the impact effect for each variable looks like. Report plots of the
      % impact effects and provide the numerical values for the cut-off points.

    ndraws=1000000;              %number of MH iterations
    nburn=500000;               %number of burn-in draws
    hmax=18;         %impulse response horizon (17 months)

    data=[GDP DEF FED M2];
    % Sign Restrictions
    yy=data;
    n=size(yy,2);
    % Get data matrices
    [X,y,T]=getXy(yy,ndet,nlags);
    yyy=y;
    xxx=X;

    bounds = 5;
    %x1=-bounds:.0001:0;    %Negative Parameters 
    %y1=0:.0001:bounds;     %Positive Parameters
    z1=-bounds:.01:bounds;
    cc=0;
    ssigma=100;
    nu=3;
    z1=-bounds:.01:bounds;  %grid for parameters where no sign is imposed a priori
    yp1=0:0.0001:0.0258;
    a_alpha = -5;
    b_alpha = 5;
    c_alpha = 0.5*(a_alpha+b_alpha);

%   Alpha1 (Uniform)
    prior_alpha1 = uniform_prior(yp1,a_alpha,b_alpha);
%   Alpha2 (Uniform)
    prior_alpha2 = uniform_prior(yp1,a_alpha,b_alpha);
%   Alpha3 (Uniform)
    prior_alpha3 = uniform_prior(yp1,a_alpha,b_alpha);
%   Beta1 (Uniform)
    prior_beta1 = uniform_prior(yp1,a_alpha,b_alpha);
%   Beta2 (Uniform)
    prior_beta2 = uniform_prior(yp1,a_alpha,b_alpha);
%   Beta3 (Uniform)
    prior_beta3 = uniform_prior(yp1,a_alpha,b_alpha);
%   Theta1 (Uniform)
    prior_theta1 = uniform_prior(yp1,a_alpha,b_alpha);
%   Theta2 (Uniform)
    prior_theta2 = uniform_prior(yp1,a_alpha,b_alpha);
%   Theta3 (Uniform)
    prior_theta3 = uniform_prior(yp1,a_alpha,b_alpha);
%   Delta1 (Uniform)
    prior_delta1 = uniform_prior(yp1,a_alpha,b_alpha);
%   Delta2 (Uniform)
    prior_delta2 = uniform_prior(yp1,a_alpha,b_alpha);     
%   Delta3 (Uniform)
    prior_delta3 = uniform_prior(yp1,a_alpha,b_alpha);

    % Set arbitrary initial values for elements of A 
    A_old=zeros(12,1);
    c=size(A_old,1);                             %number of parameters to be estimated in A

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1b: Set informative priors on lagged coefficients (B) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Compute standard deviation of each series residual via an OLS regression
% to be used in setting the prior (here: AR(24))
[s11,uhat1]=sd_prior(yy(:,1),nlags);
[s22,uhat2]=sd_prior(yy(:,2),nlags);
[s33,uhat3]=sd_prior(yy(:,3),nlags);
[s44,uhat4]=sd_prior(yy(:,4),nlags);


lambda0=10^9;    %overall confidence in prior (smaller lambda0 implies greater weight to RW)
lambda1=1;       %confidence on higher-order lags (lambda1 = 0 gives all lags equal weight)
lambda2=1;       %confidence in other-than-own lags 
lambda3=100;     %tightness of constant term 

% Specify the prior mean of the coefficients of the 3 equations of the VAR
% and their prior covariance

% EXPECTED VALUE OF REDUCED-FORM COEFFICIENTS
eta=[eye(n) zeros(n,(n*nlags+1)-n)];


% PRIOR COVARIANCE 
SS=[s11;s22;s33;s44];
M=getM(lambda0,lambda1,lambda2,lambda3,SS,nlags,n);
M1=M(:,:,1);
M2=M(:,:,2);
M3=M(:,:,3);
M4=M(:,:,4);

% Compute summary statistics of the observed data:
Syy=yyy'*yyy;
Sxx=xxx'*xxx;
Sxy=yyy'*xxx;
omega_hatT=(Syy-Sxy*inv(Sxx)*Sxy')/T;
dd3=sqrt(omega_hatT(3,3));

%uniform prior for alpha_yp/det(A) over [-1.5;0]
a_h = -1.5;
b_h = 0;
xp1=a_h:0.0001:b_h;
prior_h23 = uniform_prior(xp1,a_h,b_h);
%plot(xp1,prior_h23,'b','linewidth',3)

% Compute M_star(i) 
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

kappa=0.5;
kappastar = kappa + (T/2); % posterior kappastar

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% JIM's PROPOSAL to deal with scaling issue
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%fixed parameter values
%param=[a_alpha;b_alpha;cc;ssigma;nu;cc;ssigma;nu;a_alpha;b_alpha;cc;ssigma;nu;cc;ssigma;nu;a_alpha;b_alpha;a_alpha;b_alpha; ...
 %      cc;ssigma;nu;cc;ssigma;nu;cc;ssigma;nu;cc;ssigma;nu;T;vec(omega_hatT);kappastar;kappa];

param=[a_alpha;b_alpha;T;vec(omega_hatT);kappastar;kappa];

f_anon = @(theta_hat)post_val_KMuniform(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M4,M_star1,M_star2,M_star3,M_star4,S);

% find posterior mode
theta_zero = A_old;     %mode of prior distributions
options = optimset('LargeScale','off','MaxFunEvals',10000);
[theta_max,val_max,exitm,om,gm,HM] = fminunc(f_anon,theta_zero,options);
if min(eig(inv(HM)))>0
    PH=chol(inv(HM))';
else
    PH=eye(c);
end
PH;

%start MH algorithm with theta_hat
A_old=-theta_max;
A=[1 -A_old(1:3,1)'; -A_old(4,1) 1 -A_old(5:6,1)'; -A_old(7:8,1)' 1 -A_old(9,1); -A_old(10:12,1)' 1];

% STEP 2: Set the variance of the candidate generating density (P) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xsi=0.0002;        %tuning parameter to achieve a target acceptance rate of around 30%
W=xsi*eye(c);      %variance of RW-MH  

% STEP 3: Evaluate posterior at starting value for A:  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
omega=A*S*A';
taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);
       
% Evaluate prior p(A) at old draw
prior_a = 1/(b_alpha-a_alpha);
prior_j = 1/(b_alpha-a_alpha);
prior_k = 1/(b_alpha-a_alpha);
prior_l = 1/(b_alpha-a_alpha);
prior_e = 1/(b_alpha-a_alpha);
prior_b = 1/(b_alpha-a_alpha);
prior_c = 1/(b_alpha-a_alpha);
prior_d = 1/(b_alpha-a_alpha);
prior_f = 1/(b_alpha-a_alpha);
prior_g = 1/(b_alpha-a_alpha);
prior_h = 1/(b_alpha-a_alpha);
prior_i = 1/(b_alpha-a_alpha);

% Compute posterior value at candidate draw
log_priors=log(prior_a)+log(prior_b)+log(prior_c)+log(prior_d)+log(prior_e)+log(prior_f)+log(prior_g)+log(prior_h)+log(prior_i)+log(prior_j)+log(prior_k)+log(prior_l);

up=log_priors + T/2*log(det(A*omega_hatT*A'));
down=kappastar*log((2/T)*taustar1) ...
    +kappastar*log((2/T)*taustar2) ...
    +kappastar*log((2/T)*taustar3);
posteriorOLD=up-down;

% RW-MH algorithm 
naccept=0;
count=0;

% Store posterior distribution (after burn-in)
A_post=zeros(c+1,ndraws-nburn);
IRF=zeros(n,n,hmax,ndraws-nburn);

while count<ndraws 
      count=count+1;
      if (count/10000) == floor(count/10000)
          count
      end

    % STEP 4a: Generate draw for A from the RW candidate density %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A_new=A_old+chol(W)'*PH*randn(c,1)/sqrt(0.5*(randn(1)^2 + randn(1)^2));    % fat tails

    % STEP 4b: Evaluate posterior at new draw %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    A=[1 -A_new(1:3,1)'; -A_new(4,1) 1 -A_new(5:6,1)'; -A_new(7:8,1)' 1 -A_new(9,1); -A_new(10:12,1)' 1 ];
    H=dd3*A_new(2,1)/det(A);
    
   % if A_new(1,1)>0 && A_new(1,1)<0.0258 && A_new(2,1)<0 && A_new(3,1)<0 && ...
   %         A_new(4,1)>0 && H<0 && H>-1.5

   if A_new(1,1)>0 && A_new(1,1)<0.0258

       % Evaluate prior p(A) at new draw
        prior_a = 1/(b_alpha-a_alpha);
        prior_j = 1/(b_alpha-a_alpha);
        prior_k = 1/(b_alpha-a_alpha);
        prior_l = 1/(b_alpha-a_alpha);
        prior_e = 1/(b_alpha-a_alpha);
        prior_b = 1/(b_alpha-a_alpha);
        prior_c = 1/(b_alpha-a_alpha);
        prior_d = 1/(b_alpha-a_alpha);
        prior_f = 1/(b_alpha-a_alpha);
        prior_g = 1/(b_alpha-a_alpha);
        prior_h = 1/(b_alpha-a_alpha);
        prior_i = 1/(b_alpha-a_alpha);
       
       % Construct full matrix A     
       omega=A*S*A';
       taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
       taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
       taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
       taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);
       
       % Compute posterior value at new candidate draw
       log_priors=log(prior_a)+log(prior_b)+log(prior_c)+log(prior_d)+log(prior_e)+log(prior_f)+log(prior_g)+log(prior_h)+log(prior_i)+log(prior_j)+log(prior_k)+log(prior_l);
       up=log_priors + T/2*log(det(A*omega_hatT*A'));
       down=kappastar*log((2/T)*taustar1) ...
           +kappastar*log((2/T)*taustar2) ...
           +kappastar*log((2/T)*taustar3);
       posteriorNEW=up-down;
                
       % STEP 5: Compute acceptance probability %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       accept=min([exp(posteriorNEW-posteriorOLD);1]);
       u=rand(1);                      %draw from a uniform distribution
       if u<=accept
          A_old=A_new;                %we retain the new draw
          posteriorOLD=posteriorNEW;
          naccept=naccept+1;          %count the number of acceptances
       end
       
    end
           
    if count>nburn
           
       %Store results after burn-in          
       AA=[1 -A_old(1:3,1)'; -A_old(4,1) 1 -A_old(5:6,1)'; -A_old(7:8,1)' 1 -A_old(9,1); -A_old(10:12,1)' 1 ];
       A_post(:,count-nburn)=[A_old;dd3*A_old(2,1)/det(AA)];
       
       % STEP 7: Generate a draw for d(ii)^-1 from independent gamma %
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       omega=AA*S*AA';
       d11=inv(gamrnd(kappastar,1/gettau(kappa,omega(1,1),AA(1,:),Syy,Sxy,eta,M1,M_star1)));
       d22=inv(gamrnd(kappastar,1/gettau(kappa,omega(2,2),AA(2,:),Syy,Sxy,eta,M2,M_star2)));
       d33=inv(gamrnd(kappastar,1/gettau(kappa,omega(3,3),AA(3,:),Syy,Sxy,eta,M3,M_star3)));
       d44=inv(gamrnd(kappastar,1/gettau(kappa,omega(4,4),AA(4,:),Syy,Sxy,eta,M4,M_star4)));
       DD=diag([d11;d22;d33;d44]);
       
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
jd=3;
[agi,bgi]=hist((1./A_post(jd,:)),10000000);
deltai=bgi(1,2)-bgi(1,1);
post_ii=agi./((ndraws-nburn)*deltai);

% Compute acceptance ratio of RW-MH algorithm
acceptance_ratio=naccept/ndraws;
disp(['Acceptance ratio:' num2str(acceptance_ratio)])


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plot impulse responses for one-standard deviation shocks
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%seriesnames=["GDP","GDP Deflactor","M2","FFR"];

% alph=0.16;
% index=[alph*(ndraws-nburn) (1-alph)*(ndraws-nburn)];    %implies 95% coverage of the entire distribution
% HO=(0:1:hmax-1)';                                       %impulse response horizon
% figure(6)
% subplot(4,4,1)
% x1=-sort(cumsum(squeeze(IRF(1,1,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('GDP','fontsize',14)
% title('Demand shock','fontsize',14)
% 
% subplot(4,4,5)
% x1=sort(cumsum(squeeze(IRF(1,2,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('GDP','fontsize',14)
% title('Supply shock','fontsize',14)
% 
% subplot(4,4,9)
% x1=sort(cumsum(squeeze(IRF(1,3,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('GDP','fontsize',14)
% title('MP shock','fontsize',14)
% %xlabel('Months','fontsize',12)
% 
% subplot(4,4,13)
% x1=sort(cumsum(squeeze(IRF(1,4,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('GDP','fontsize',14)
% title('M2 shock','fontsize',14)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(4,4,2)
% x1=-sort(cumsum(squeeze(IRF(2,1,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('Inflation','fontsize',14)
% title('Demand shock','fontsize',14)
% 
% subplot(4,4,6)
% x1=sort(cumsum(squeeze(IRF(2,2,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('Inflation','fontsize',14)
% title('Supply shock','fontsize',14)
% 
% subplot(4,4,10)
% x1=sort(cumsum(squeeze(IRF(2,3,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('Inflation','fontsize',14)
% title('MP shock','fontsize',14)
% %xlabel('Months','fontsize',12)
% 
% subplot(4,4,14)
% x1=sort(cumsum(squeeze(IRF(2,4,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis([0 hmax-1 -2.5 1.5])
% ylabel('Inflation','fontsize',14)
% title('M2 shock','fontsize',14)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(4,4,3)
% x1=-sort(cumsum(squeeze(IRF(3,1,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('FFR','fontsize',14)
% title('Demand shock','fontsize',14)
% 
% subplot(4,4,7)
% x1=sort(cumsum(squeeze(IRF(3,2,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('FFR','fontsize',14)
% title('Supply shock','fontsize',14)
% 
% subplot(4,4,11)
% x1=sort(cumsum(squeeze(IRF(3,3,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('FFR','fontsize',14)
% title('MP shock','fontsize',14)
% %xlabel('Months','fontsize',12)
% 
% subplot(4,4,15)
% x1=sort(cumsum(squeeze(IRF(3,4,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('FFR','fontsize',14)
% title('M2 shock','fontsize',14)
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% subplot(4,4,4)
% x1=-sort(cumsum(squeeze(IRF(4,1,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('M2','fontsize',14)
% title('Demand shock','fontsize',14)
% 
% subplot(4,4,8)
% x1=sort(cumsum(squeeze(IRF(4,2,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('M2','fontsize',14)
% title('Supply shock','fontsize',14)
% 
% subplot(4,4,12)
% x1=sort(cumsum(squeeze(IRF(4,3,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('M2','fontsize',14)
% title('MP shock','fontsize',14)
% %xlabel('Months','fontsize',12)
% 
% subplot(4,4,16)
% x1=sort(cumsum(squeeze(IRF(4,4,:,:)),1),2);
% temp1=[(median(x1,2)) x1(:,index(1)) x1(:,index(2))];
% plotx1(temp1,HO); box on; plot(HO,zeros(hmax,1),'k:')
% axis tight
% ylabel('M2','fontsize',14)
% title('M2 shock','fontsize',14)


%%IMPACT EFFECT USING SIGN RESTRICTIONS IN A BAYESIAN APPROACH

figure(7)

subplot(4,4,1)
hist(squeeze(IRF(1,1,1,:)),500)
axis tight
ylabel('GDP','fontsize',14)
title('Demand shock','fontsize',14)

subplot(4,4,5)
hist(squeeze(IRF(1,2,1,:)),500)
axis tight
ylabel('GDP','fontsize',14)
title('Supply shock','fontsize',14)

subplot(4,4,9)
hist(squeeze(IRF(1,3,1,:)),500)
axis tight
ylabel('GDP','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,13)
hist(squeeze(IRF(1,4,1,:)),500)
axis tight
ylabel('GDP','fontsize',14)
title('M2 shock','fontsize',14)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,4,2)
hist(squeeze(IRF(2,1,1,:)),500)
axis tight
ylabel('Inflation','fontsize',14)
title('Demand shock','fontsize',14)

subplot(4,4,6)
hist(squeeze(IRF(2,2,1,:)),500)
axis tight
ylabel('Inflation','fontsize',14)
title('Supply shock','fontsize',14)

subplot(4,4,10)
hist(squeeze(IRF(2,3,1,:)),500)
axis tight
ylabel('Inflation','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,14)
hist(squeeze(IRF(2,4,1,:)),500)
axis tight
ylabel('Inflation','fontsize',14)
title('M2 shock','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,4,3)
hist(squeeze(IRF(3,1,1,:)),500)
axis tight
ylabel('FFR','fontsize',14)
title('Demand shock','fontsize',14)

subplot(4,4,7)
hist(squeeze(IRF(3,2,1,:)),500)
axis tight
ylabel('FFR','fontsize',14)
title('Supply shock','fontsize',14)

subplot(4,4,11)
hist(squeeze(IRF(3,3,1,:)),500)
axis tight
ylabel('FFR','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,15)
hist(squeeze(IRF(3,4,1,:)),500)
axis tight
ylabel('FFR','fontsize',14)
title('M2 shock','fontsize',14)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
subplot(4,4,4)
hist(squeeze(IRF(4,1,1,:)),500)
axis tight
ylabel('M2','fontsize',14)
title('Demand shock','fontsize',14)

subplot(4,4,8)
hist(squeeze(IRF(4,2,1,:)),500)
axis tight
ylabel('M2','fontsize',14)
title('Supply shock','fontsize',14)

subplot(4,4,12)
hist(squeeze(IRF(4,3,1,:)),500)
axis tight
ylabel('M2','fontsize',14)
title('MP shock','fontsize',14)
%xlabel('Months','fontsize',12)

subplot(4,4,16)
hist(squeeze(IRF(4,4,1,:)),500)
axis tight
ylabel('M2','fontsize',14)
title('M2 shock','fontsize',14)

%%Finding the numerical values for the cut-off points%%%
%Sign-restriction:
max_fedrate_sign = [ max(squeeze(IRF(1,1,1,:))) max(squeeze(IRF(1,2,1,:))) max(squeeze(IRF(1,3,1,:))) max(squeeze(IRF(1,4,1,:)))];
min_fedrate_sign = [ min(squeeze(IRF(1,1,1,:))) min(squeeze(IRF(1,2,1,:))) min(squeeze(IRF(1,3,1,:))) min(squeeze(IRF(1,4,1,:)))];
max_M2_sign = [ max(squeeze(IRF(2,1,1,:))) max(squeeze(IRF(2,2,1,:))) max(squeeze(IRF(2,3,1,:))) max(squeeze(IRF(2,4,1,:)))];
min_M2_sign = [ min(squeeze(IRF(2,1,1,:))) min(squeeze(IRF(2,2,1,:))) min(squeeze(IRF(2,3,1,:))) min(squeeze(IRF(2,4,1,:)))];
max_inflation_sign = [ max(squeeze(IRF(3,1,1,:))) max(squeeze(IRF(3,2,1,:))) max(squeeze(IRF(3,3,1,:))) max(squeeze(IRF(3,4,1,:)))];
min_inflation_sign = [ min(squeeze(IRF(3,1,1,:))) min(squeeze(IRF(3,2,1,:))) min(squeeze(IRF(3,3,1,:))) min(squeeze(IRF(3,4,1,:)))];
max_GDP_sign = [ max(squeeze(IRF(4,1,1,:))) max(squeeze(IRF(4,2,1,:))) max(squeeze(IRF(4,3,1,:))) max(squeeze(IRF(4,4,1,:)))];
min_GDP_sign = [ min(squeeze(IRF(4,1,1,:))) min(squeeze(IRF(4,2,1,:))) min(squeeze(IRF(4,3,1,:))) min(squeeze(IRF(4,4,1,:)))];

% Table
Puntos_corte_max = [max_GDP_sign; max_inflation_sign; max_M2_sign; max_fedrate_sign];
Puntos_corte_min = [min_GDP_sign; min_inflation_sign; min_M2_sign; min_fedrate_sign];

names = {'Demand shock';'Supply shock';'MP shock';'M2 shock'};
names = convertCharsToStrings(names);
series = {'Demand';'Inflation';'FFR';'M2'};
series=cellstr(series);


Puntos_corte_max = array2table(Puntos_corte_max); Puntos_corte_max.Properties.VariableNames = names;
Puntos_corte_max.Properties.RowNames = series;

Puntos_corte_min = array2table(Puntos_corte_min); Puntos_corte_min.Properties.VariableNames = names;
Puntos_corte_min.Properties.RowNames = series;

%Saving the Tables
%writetable(Puntos_corte_max, 'C:\Users\danil\OneDrive - Universidad del Pacífico\BCRP\A Bayesian Approach to Identification of Structural VAR Models\Trabajo\Resultados\Puntos_corte.xlsx','sheet','MAX','WriteRowNames',true)
%writetable(Puntos_corte_min, 'C:\Users\danil\OneDrive - Universidad del Pacífico\BCRP\A Bayesian Approach to Identification of Structural VAR Models\Trabajo\Resultados\Puntos_corte.xlsx','sheet','MIN','WriteRowNames',true)
