function [Vlambda,gammahat] = GMM_variance(X,ehat,z)
% inputs
%     X = (T x k) matrix of predictor variables for VAR
%     ehat = (T x n) matrix of OLS residuals
%     z = (Ts x 1) vector of scalar instruments beginning with obs (T - Ts + 1) of ehat
% outputs
%     gammahat = (n x 1) vector of coefficients from regressions of ehat on z over the (Ts x 1) subsample
%     Vlambda = T times asymptotic variance of (pi,gammahat)

% ===========================================================
% calculate gammahat

T = size(ehat,1);
n = size(ehat,2);
k = size(X,2);
Ts = size(z,1);

if T > Ts
    zeroz = zeros(T-Ts,1);
    z = [zeroz; z];
end

gammahat = (ehat'*z)/(z'*z);

zresids = (ehat - (z*gammahat')).*(z*ones(1,n));
    % zresids is (T x n)

% ===========================================================
% calculate Qhat

nk = n*k;
Qhat = zeros(nk+n,nk+n);
for t =1:T
   qt = [kron(ehat(t,:),X(t,:)) zresids(t,:)];
   Qhat = Qhat + qt'*qt;
end

Qhat(1:nk,1:nk) = Qhat(1:nk,1:nk)/T;
Qhat(nk+1:end,:) = Qhat(nk+1:end,:)/Ts;
Qhat(1:nk,nk+1:end) = Qhat(1:nk,nk+1:end)/Ts;

% ===========================================================
% calculate inv(D')

Qx = inv((X'*X)/T);

sz_inv = Ts/sum(z.^2);
Qzx = z'*X/Ts;
invDprime = [kron(eye(n),-Qx) zeros(nk,n); ...
             kron(eye(n),sz_inv*Qzx*Qx) (-eye(n)*sz_inv)];

% ===========================================================
% calculate Vlambda

Shat = Qhat;
Shat(nk+1:end,nk+1:end) = Shat(nk+1:end,nk+1:end)*T/Ts;
Vlambda = invDprime*Shat*invDprime';

end