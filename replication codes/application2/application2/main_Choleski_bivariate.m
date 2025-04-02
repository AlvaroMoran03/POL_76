%%% Traditional sign restriction algorithm
%%% Bivariate oil market VAR model without sign restrictions
%%% By Christiane Baumeister, October 2017

clear; clc;
load oildata.txt
% oil production, real economic activity, real price of oil

quantity=oildata(:,1);
price=oildata(:,3);
data=[quantity price];
                   
const=ones(size(data,1),1);   %deterministic term: constant
ndet=1;                       %number of deterministic variables
data=[const data];

nsteps=17;                %horizon for impulse responses                          
   
%model setup
[sp,nvars]=size(data);   %sp is the total number of observations
nvars=nvars-ndet;        %the '-ndet' takes out the counting of deterministic variables
nlags=24;                %number of lags

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


%point estimation
P=chol(vmat)';     % Choleski decomposition
          
%compute structural impulse responses:
irf=impulse(beta',P,ndet,nlags,nsteps);        
   % se necesita la forma reducida para poder computar el impulso respuesta

% Impulse responses of one-standard-deviation shocks
% oil supply shock
irsup1=squeeze(irf(1,1,:))';   %response of variable 1
irsup2=squeeze(irf(2,1,:))';   %response of variable 2
% oil demand shock
irdem1=squeeze(irf(1,2,:))';   %response of variable 1
irdem2=squeeze(irf(2,2,:))';   %response of variable 2

% cumulate production growth rates to levels
irsup1=cumsum(irsup1,2);
irdem1=cumsum(irdem1,2);


% graphs of structural IRFs

h=(0:1:nsteps)';

figure(1)
subplot(221)
plot(h,-irsup1,'b','linewidth',2); grid;
axis([0 16 -2 1])
title('Oil production')
ylabel('Oil supply shock')
subplot(223)
plot(h,-irsup2,'b','linewidth',2); grid;
axis([0 16 -5 10])
title('Real oil price')
ylabel('Oil supply shock')

subplot(222)
plot(h,irdem1,'b','linewidth',2); grid;
axis([0 16 -2 1])
title('Oil production')
ylabel('Oil demand shock')
subplot(224)
plot(h,irdem2,'b','linewidth',2); grid; 
axis([0 16 -5 10])
title('Real oil price')
ylabel('Oil demand shock')


%compute structural parameters
D=eye(2);
D(1,1)=P(1,1)^2;
D(2,2)=P(2,2)^2;
A=D^.5*inv(P)     %matrix of structural parameters

%Which structural parameter is this (economic interpretation)?
z=-A(2,2)/A(2,1);

