%%% Application 2: Bivariate Labor Market VAR 

%%% Identification of 2 shocks with sign restrictions: 
%%% demand shock (+,+) and supply shock (+,-)

%%% By Christiane Baumeister and James D. Hamilton (Sept 2013)

clear; clc;

seednumber=140778;
rand('seed',seednumber);
randn('seed',seednumber);

load labor_data               %frequency: quarterly
time=(1947.25:0.25:2014.5)';  %1947Q1 to 2014Q2
d=[time data];
data=data(84:end,:);          %1967Q4: 84
time=time(84:end,1);
wage=100*log(data(:,1));      %real compensation per hour (from FRED)
empl=100*log(data(:,2));      %total employment

clear data
yy=[lagn(wage,1) lagn(empl,1)]; %quarterly growth rates
time=time(2:end,1);

size(yy)
nlags=8;                      %number of lags

data=yy;
const=ones(size(data,1),1);   %deterministic term: constant
trend=(1:1:size(data,1))';    %time trend
ndet=1;                       %number of deterministic variables

if ndet==1
    data=[const data];
elseif ndet==2
    data=[const trend data];
end

nsteps=1;                 %horizon for impulse responses                                             
ndraws=20000;             %number of solutions that match the sign restrictions
bmed=0.5;                 %impulse response percentile (median)
bupp=0.84;                %upper error band (percentile)
blow=0.16;                %lower error band (percentile)
cors=1;                   %first period sign restrictions are imposed
h=1;                      %maximum length for which the sign restrictions are imposed    
    
[sp,nvars]=size(data);   %sp is the total number of observations
nvars=nvars-ndet;        %takes out the counting of deterministic variables
                 

ess=sp-nlags;       %effective sample size after taking lags into account
sb=nlags+1;         %sample beginning
sl=sp;              %last period of sample
ncoe=nvars*nlags;   %number of coefficients without deterministic variables

%construct X for Y=X*B+U

x=data(:,ndet+1:end);

X=zeros(ess,ncoe);
for k=1:nlags
    X(:,nvars*(k-1)+1:nvars*k) = x(sb-k:sl-k,:);  %first lag of all variables, second lag of all variables etc.
end

if ndet==1
    X=[ones(ess,1) X];   %without trend
elseif ndet==2
    X=[ones(ess,1) trend(1:ess,1) X];
end

%estimation

y=x(sb:sl,:);
xtx=X'*X;
xty=X'*y;
%OLS estimation of reduced-form coefficients
%Bh=X\(X'\xty);     %X\X' is a faster alternative for inv(X'X)
beta=inv(X'*X)*X'*y;
e=y-X*beta;
%variance-covariance matrix 
vmat=y'*(eye(ess)-X*inv(X'*X)*X')*y/ess;
%this way of calculating sigma yields the same result as the built-in in RATS
%alternative way of calculation taking correction for degrees of freedom
%into account:
%vmat=y'*(eye(ess)-X*inv(xtx)*X')*y/(ess-ncoe-1);
% [A,SIGMA,U,V,XX]=olsvarc(yy,nlags);
% vmat=SIGMA;

% Bayesian estimation
% kfset=inv(X'*X);
% sxx=chol(kfset)';
% svtr=chol(vmat);
% betaols=beta;
% ncoefs=size(sxx,1);


%Monte Carlo simulation
tel=0;
count=0;

irsup1=zeros(ndraws,nsteps);
irsup2=zeros(ndraws,nsteps);

irdem1=zeros(ndraws,nsteps);
irdem2=zeros(ndraws,nsteps);

disp('Starting Monte Carlo Simulations')
while tel<ndraws
    count=count+1;
    
%     BAYESIAN: draws from posterior
%     temp01=size(y,1)*inv(wishrnd(eye(nvars),size(y,1)));   %inverse Wishart
%     sigmad=svtr'*temp01*svtr;
%     swish=chol(sigmad)';
%     ranc=randn(ncoefs,nvars);
%     betau=sxx*ranc*swish';                                 %normal 
%     betadraw=betaols+betau; %contains the posterior draws for the coefficients in the following form: y(t)=c+y(t-1)+y(t-2)+y(t-3)...
                              %where y contains all the variables
                              %so, first column: constant, coefficient first
                              %lag first variable, first lag second
                              %variable, first lag third variable...,
                              %second lag first variable, second lag second
                              %variable...,
    
        
    %generate rotation matrix by QR decomposition of draw from
    %multivariate normal (Rubio-Ramirez, Waggoner and Zha, REStud 2010)    
    K=normrnd(0,1,nvars,nvars);   %draws an nvars x nvars matrix from a standard normal distribution
    [Q,R]=qr(K);                  %QR decomposition: Q is an orthogonal matrix
    
    for i=1:nvars                 
        if R(i,i)<0
            Q(:,i)=-Q(:,i);
        end
    end
    
    %candidate impact matrix:
    CI=chol(vmat)'*Q;     %for point estimation
    %CI=chol(sigmad)'*Q;  %for Bayesian estimation 
    
    %compute candidate impulse responses:
    irf=impulse(beta',CI,ndet,nlags,nsteps);       %for point estimation
    %irf=impulse(betadraw',CI,ndet,nlags,nsteps);  %for Bayesian estimation
   
    %ir=cumsum(irf,3);   %gives the accumulated impulse responses (for the
                         %case of first difference estimation)
    ncorr=nvars*nvars;            %number of sign correlations
    co=zeros(cors+h,ncorr);
    
    %normalization so that you check for positive and negative shocks
    for sh=1:nvars
        for k=1:cors+h
            co(k,1+nvars*(sh-1))=irf(1,sh,k)/irf(1,sh,1); 
            co(k,2+nvars*(sh-1))=irf(2,sh,k)/irf(1,sh,1);
        end
    end
    
    s1is=0;
    s2is=0;
   
    %%%Sign restrictions only on impact
        
        % Aggregate supply shock (label 100)
        if co(1,1)>=0 && co(1,2)<=0
            s1is=100; %Le asigno etiquetas
        end
        if co(1,3)>=0 && co(1,4)<=0
            s2is=100;
        end
                
        % Aggregate demand shock (label 10)
        if co(1,1)>=0 && co(1,2)>=0 
            s1is=10;
        end
        if co(1,3)>=0 && co(1,4)>=0 
            s2is=10;
        end
                
        
    identot=s1is+s2is;
    clear co
    
    if identot==110 % Es un supply shock
        tel=tel+1
       
        %shock 1
  
        if s1is==100  %then it is a supply shock
               irsup1(tel,1:nsteps)=irf(1,1,1:nsteps);
               irsup2(tel,1:nsteps)=irf(2,1,1:nsteps);
        end
        
        if s1is==10   %then it is a demand shock
               irdem1(tel,1:nsteps)=irf(1,1,1:nsteps);
               irdem2(tel,1:nsteps)=irf(2,1,1:nsteps);
        end
        
         
       %shock 2
  
        if s2is==100
                irsup1(tel,1:nsteps)=irf(1,2,1:nsteps);
                irsup2(tel,1:nsteps)=irf(2,2,1:nsteps);
        end
        
        if s2is==10
               irdem1(tel,1:nsteps)=irf(1,2,1:nsteps);
               irdem2(tel,1:nsteps)=irf(2,2,1:nsteps);
        end
        
       
    end
    
end

disp(['Number of full identifications: ' num2str(tel)])
disp(['Total number of draws and rotations: ' num2str(count)])

%normalization of the signs to compute averages, min and max

tester2=irdem1(:,1);
tester3=irsup1(:,1);   

vts2=zeros(ndraws,1);
vts3=zeros(ndraws,1);

for jx=1:tel
     
    if tester2(jx,1)>0
        vts2(jx,1)=1;
    else vts2(jx,1)=-1;
    end
        
    if tester3(jx,1)>0
        vts3(jx,1)=1;
    else vts3(jx,1)=-1;
    end
end

%impulse responses (adjust for positive and negative shocks)
% Ajusta si es positiva o negativa
for jy=1:nsteps
    irsup1(:,jy)=vts3.*irsup1(:,jy);
    irsup2(:,jy)=vts3.*irsup2(:,jy);
end

for jy=1:nsteps
    irdem1(:,jy)=vts2.*irdem1(:,jy);
    irdem2(:,jy)=vts2.*irdem2(:,jy);
end

%normalization
% Hago la matriz del 72 pag 
irsup1_norm=irsup1(:,1)./irsup1(:,1);
irsup2_norm=irsup2(:,1)./irsup1(:,1);

irdem1_norm=irdem1(:,1)./irdem1(:,1);
irdem2_norm=irdem2(:,1)./irdem1(:,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
g=irdem2_norm;    %supply elasticity
h=irsup2_norm;    %demand elasticity

%%% Plots for g
nbin=500;
[ag,bg]=hist(g,nbin);
delta=bg(1,2)-bg(1,1);
%delta=(max(h)-min(h))/nbin;   %equivalent
aag=ag./(ndraws*delta);

vmatc=chol(vmat)';
p11=vmatc(1,1);
p21=vmatc(2,1);
p22=vmatc(2,2);
gamm=p22/p11;   %scale
x0=p21/p11;     %location

gH=vmat(2,2)/vmat(1,2)
gL=vmat(1,2)/vmat(1,1)
dd=(gL:0.0001:gH)';

FH=0.5+pi\atan((gH-x0)/gamm);
FL=0.5+pi\atan((gL-x0)/gamm);
fxg=(1/pi.*(gamm./((dd-x0).^2+gamm^2)))/(FH-FL);

% f=@(dd) (1/pi.*(gamm./((dd-x0).^2+gamm^2)))/(FH-FL);
% Q=integral(f,gL,gH)

%%% Plots for h 
ss=(-1000000:0.1:0)';
F0=0.5+pi\atan((0-x0)/gamm);
fxh=(1/pi.*(gamm./((ss-x0).^2+gamm^2)))/F0;
[ah,bh]=hist(h,ndraws*10);
delta=bh(1,2)-bh(1,1);
aah=ah./(ndraws*delta);

fh=@(ss) (1/pi.*(gamm./((ss-x0).^2+gamm^2)))/F0;
Q=integral(fh,-inf,0)


figure(1)
subplot(2,2,1)
bar(bh,aah), hold on, plot(ss,fxh,'r','linewidth',2)
axis([-5 0 0 2])
title('Demand elasticity','fontsize',12)

subplot(2,2,3)
bar(bg,aag), hold on, plot(dd,fxg,'r','linewidth',2), ...
    hold on, line([gL gL], [0 2]), ...
    hold on, line([gH gH], [0 2]), ...
    hold on, plot(-ss,fxh,'g','linewidth',2),
axis([0 5 0 2])
text(0.075,1.85,'\ith_{\itL}','fontsize',10,'color','r')
text(4.1,1.85,'\ith_{\itH}','fontsize',10,'color','r')
title('Supply elasticity','fontsize',12)

