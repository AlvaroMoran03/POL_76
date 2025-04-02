% This Matlab code calculates the structural coefficients of the Taylor
% rule using the approach proposed in Baumeister and Hamilton (ET 2024)

clear;
clc;

% load data
var3_data = csvread('var3_data.csv');
   % var3_data is (275 x 4) matrix
   % row 1 = 1954:Q3, row 275 is 2023:Q1
   % col 1 = date
   % col 2:  level of real GDP (2012 dollars)
   % col 3:  level of implicit price deflator on personal consumption expenditures (2009 = 100)
   % col 4:  average value of effective fed funds rate or Wu-Xia shadow rate over the quarter
   n=size(var3_data,2)-1;

% convert GDP and price level to annual growth rates
  gdp_growth = 400*(log(var3_data(2:end,2)) - log(var3_data(1:end-1,2)));
  inflation = 400*(log(var3_data(2:end,3)) - log(var3_data(1:end-1,3)));
  Y = [gdp_growth inflation var3_data(2:end,4)];
  % Y is (274 x 3) matrix, row 1 = 1954:Q4 row 261 = 2019:Q4

  % order interest rate first in VAR
  Y = [Y(:,3) Y(:,1:2)];

  % end data with 2019:Q4
  Y = Y(1:261,:);

    instdata = csvread('Bauer_Swanson2.csv');
    % col 1 = date, col 2 = sum over three months of quarter of Bauer-Swanson instrument
    % row 1 = 1988:Q1, row 128 = 2019:Q4
    z = instdata(:,2);

% ==================================================================================
% estimate VAR

nlags = 4;   % number of lags in VAR
T = size(Y,1)-nlags;  % T is the number of observations used for estimation
ndet=1;       % constant included in VAR
[X,y]=getXy(Y,ndet,nlags);
   % y is (T x n) = (257 x 3)
   % X is (T x k) for k = 1 + m*n = 13

Pihat=(X'*X)\X'*y; % Pihat is (k x 3)
ehat=y-X*Pihat;    % ehat is (T x n)
vmat=1/size(ehat,1)*(ehat'*ehat); % variance-covariance matrix

k = size(X,2);
Ts = size(z,1);

if T > Ts
    zeroz = zeros(T-Ts,1);
    z = [zeroz; z];
end

% impact effect of instrument on interest rate, real GDP, inflation
gammahat = (ehat'*z)/(z'*z); 
% estimate structural parameters in the Taylor rule (not normalized)
vhat = inv(vmat)*gammahat;
% alternative way:
vhat_alt = (z'*ehat)/(ehat'*ehat); % OLS regression of z on ehat
'policy rule coefficients'
'not normalized'
vhat'

vhat = vhat/vhat(1);
vhat_alt = vhat_alt/vhat_alt(1);
'normalized'
vhat'
vhat_alt

% ==================================================================================
% Calculate covariance between residuals and instruments and asymptotic variance matrix
% [Vlambda,gammahat] = GMM_variance(X,ehat,z);

% 'coeffs std errors'
% pihat = reshape(Pihat,size(Pihat,1)*size(Pihat,2),1);
% [[pihat; gammahat] sqrt(diag(Vlambda/T))]
% 
% lambdahat = [pihat; gammahat];
% Vlambdahat = Vlambda/T;
% save("gmm_results.mat","lambdahat","Vlambdahat");

% 'effects of policy'
% gammahat/gammahat(1)

% =======================================================================
% Calculate Anderson-Rubin confidence intervals
% Vgamma = Vlambda(end-n+1:end,end-n+1:end)/T;
% alpha_size = 0.10;
% 
% 'size for AR confidence interval'
% alpha_size
% 
% % interval for output effect
% alpha_1z = gammahat(1);
% alpha_iz = gammahat(2);
% sigma_ii = Vgamma(2,2);
% sigma_i1 = Vgamma(1,2);
% sigma_11 = Vgamma(1,1);
% [hupper,hlower] = AR_interval(alpha_iz,alpha_1z,sigma_ii,sigma_i1,...
%     sigma_11,alpha_size);
% 'AR interval for effect on output'
% [hupper hlower]
% 
% % interval for inflation effect
% alpha_1z = gammahat(1);
% alpha_iz = gammahat(3);
% sigma_ii = Vgamma(3,3);
% sigma_i1 = Vgamma(1,3);
% sigma_11 = Vgamma(1,1);
% [hupper,hlower] = AR_interval(alpha_iz,alpha_1z,sigma_ii,sigma_i1,...
%     sigma_11,alpha_size);
% 
% 'AR interval for effect on inflation'
% [hupper hlower]

