% Estimate univariate linear regression model with Gibbs sampling

clear;
clc;

% Generate data
% Model
% y_{t} = BETA * x_{t} + v_{t} 
% v_{t} ~ N(0,SIG2)

BETA    = 0.9;
SIG2    = 0.006;
k       = 1;   % # RHS variables
t       = 1000; % sample length
e1      = randn(t,1)*sqrt(SIG2);  %generate errors which are N(0,SIG2)
y       = zeros(t,1);
x       = randn(t,1);   %explanatory variable: i.i.d. draws from N(0,1)

%compute the dependent variable
for j=1:t
    y(j) = x(j,:)*BETA + e1(j);
end

% OLS estimation (for comparison)

bols    = inv(x'*x)*(x'*y);
sigols  = ((y-x*bols)'*(y-x*bols))/(t-k);

% BAYESIAN ESTIMATION
% Set priors for distribution of beta and sigma
% Assume    beta0 ~ N(b0,P0)
%           1/sig20 ~ G(t0,R0) where G is the Gamma distribution

b0 = 0.5;     % mean of prior distribution for beta
P0 = 10;      % variance of prior distribution for beta

% plot prior for beta
z1 = 0:0.01:1;
pn = normpdf(z1,b0,P0);
%plot(z1,pn)

t0 = 0;   % prior shape parameter for sigma
R0 = 0;   % prior scale parameter for sigma: zero implies a "flat" prior 
          % i.e. no prior information

% Set parameters for Gibbs sampler
nburn   = 90000;   % burn-in period: draws are discarded
nkeep   = 10000;   % numbers of draws on which inference will be based
n       = nkeep+nburn;  % total number of Gibbs sampling iterations

sig2    = 1;    % starting value for Gibbs sampler

Beta    = zeros(nkeep,k);
Sig2    = zeros(nkeep,1);

for iter=1:n
    
    % generate a draw for beta conditional on sig2
    P1  = inv( inv(P0) + (1/sig2)*(x'*x) );
    b1  = P1*(inv(P0)*b0 + (1/sig2)*(x'*y));
    % reject unstable draws
%     chk = -1;
%     while chk<1
    beta = b1 + chol(P1)*randn(k,1);
%         if beta<1 
%            chk=1; 
%         end
%     end
    
    % generate a draw for sig2 conditional on beta
    t1      = t0+t;
    R1      = R0+(y-x*beta)'*(y-x*beta);
    shpe    = t1;
    scle    = 1/(R1);  
    sig2inv = gamrnd(shpe,scle,1,1);
    sig2    = 1/sig2inv;
    
%     z = randn(t1,1);
%     a = z'*z;
%     sig2 = R1/a;
    
    if iter>nburn
        Beta(iter-nburn,:)  = beta;
        Sig2(iter-nburn,:)  = sig2;
    end
end

figure(1)
dim1 = 2;
dim2 = 2;
nbins = 500;

subplot(dim1,dim2,1)
hist(Beta,nbins)
title('\beta')
v=axis;
axis([0.4 1 v(3) v(4)])
line(BETA*ones(1,2),[v(3) v(4)],'Linewidth',1,'Color','r')
line(b0*ones(1,2),[v(3) v(4)],'Linewidth',1,'Color','k')
legend('estimated','true','prior')

subplot(dim1,dim2,3)
hist(Beta,nbins)
title('\beta')
v=axis;
line(BETA*ones(1,2),[v(3) v(4)],'Linewidth',1,'Color','r')
line(mean(Beta)*ones(1,2),[v(3) v(4)],'Linewidth',1,'Color','g')
legend('estimated','true','posterior mean')

subplot(dim1,dim2,4)
hist(Sig2,nbins)
title('\sigma^{2}')
v = axis;
line(SIG2*ones(1,2),[v(3) v(4)],'Linewidth',1,'Color','r')
line(mean(Sig2)*ones(1,2),[v(3) v(4)],'Linewidth',1,'Color','g')

