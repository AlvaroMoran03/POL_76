%%% Traditional sign restriction algorithm
%%% 6-variable VAR model without sign restrictions
%%% By Christiane Baumeister, October 2017

clear; clc;
load data6     %1962.Q1-2008.Q4
% 10-year yield, fed funds rate, real GDP, GDP deflator, commodity price
% index, investment

data=[data(2:end,1:2) lagn(100*log(data(:,3:end)),1)];
                   
const=ones(size(data,1),1);   %deterministic term: constant
ndet=1;                       %number of deterministic variables
data=[const data];

nsteps=1;                 %horizon for impulse responses                          
ndraws=20000;             %number of solutions that match the sign restrictions
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
%vmat=y'*(eye(ess)-X*inv(X'*X)*X')*y/ess;
%this way of calculating covariance matrix yields the same result as in RATS
%alternative way of calculation taking correction for degrees of freedom
%into account:
%sigma=y'*(eye(ess)-X*inv(xtx)*X'*y/(ess-ncoe)


% BAYESIAN setup
% kfset=inv(X'*X);
% sxx=chol(kfset)';
% svtr=chol(vmat);
% betaols=beta;
% ncoefs=size(sxx,1);


%Monte Carlo simulation

seednumber=372398;
rand('seed',seednumber);
randn('seed',seednumber);

tel=0;
count=0;

disp('Starting Monte Carlo Simulations')
while count<ndraws
    count=count+1
    
%   BAYESIAN: draws from posterior
%     temp01=size(y,1)*inv(wishrnd(eye(nvars),size(y,1)));   
%     sigmad=svtr'*temp01*svtr;
%     swish=chol(sigmad)';
%     ranc=randn(ncoefs,nvars);
%     betau=sxx*ranc*swish';
%     betadraw=betaols+betau; %contains the posterior draws for the coefficients in the following form: y(t)=c+y(t-1)+y(t-2)+y(t-3)...
                              %where y contains all the variables
                              %so, first column: constant, coefficient first
                              %lag first variable, first lag second
                              %variable, first lag third variable...,
                              %second lag first variable, second lag second
                              %variable...,
    
    %rotation matrix
    K=normrnd(0,1,nvars,nvars);   %draws an nvars x nvars matrix from a standard normal distribution
    
    [Q,R]=qr(K);                  %qr decomposition: Q is a unitary matrix
                                  %alternatively to qr decomp, you could
    for i=1:nvars                 %just calculate Q=orth(K) which delivers orthonormal matrix
        if R(i,i)<0
            Q(:,i)=-Q(:,i);
        end
    end
    
     
%     E=eig(sigmad);            %eigenvalues in a column vector
%     [eigdec,D]=eig(sigmad);   
%     
%     for je=1:nvars
%         eigdecscaled(:,je)=eigdec(:,je)*sqrt(E(je,1));   %causes each eigenvector (each column of eigdec) 
%                                                          %to have length equal to its eigenvalue
%     end
    %candidate impact matrix:
    %a0=eigdecscaled*Q;
    a0=chol(vmat)'*Q;     %for point estimation
    %a0=chol(sigmad)'*Q;  %for Bayesian estimation 
          
    %compute candidate impulse responses:
    irf(:,:,:,count)=impulse(beta',a0,ndet,nlags,nsteps);        %for point estimation
    %irf(:,:,:,count)=impulse(betadraw',a0,ndet,nlags,nsteps);   %for Bayesian estimation                        
         
end

disp(['Number of full identifications: ' num2str(tel)])
disp(['Total number of draws and rotations: ' num2str(count)])
%
%
% Impact effect of one-standard-deviation shocks (histogram)
figure(1)
subplot(221)
hist(squeeze(irf(1,1,1,:)),500)
title('10-year bond yield')
subplot(222)
hist(squeeze(irf(6,1,1,:)),500)
title('Investment')

vmat