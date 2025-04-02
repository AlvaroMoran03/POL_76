%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Cleaning the workspace
clear;
clc;

% Loading the CSV file
data = readtable('GDPC1.csv');

% Checking
% disp(head(data));

% Converting date column to datetime format
data.observation_date = datetime(data.observation_date, 'InputFormat', 'yyyy-MM-dd');

% Checking
% disp(head(data));

% Defining the date range
start_date = datetime('1948-01-01');
end_date   = datetime('2019-12-31');

% Filtering the data
filtered_data = data(data.observation_date >= start_date & data.observation_date <= end_date, :);

% Verifying the filtered data
% disp(head(filtered_data));
% disp(tail(filtered_data));

% Filtered data
df   = filtered_data;

% Calculating the annual growth rate of real GDP
GDP  = df.GDPC1;
LGDP = log(GDP);
GDPR = 100 * (LGDP(5:end) - LGDP(1:end-4)); % Annual growth rate in percent

% Dates
dates = df.observation_date(5:end);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B.1. Plot the growth rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1);
plot(dates, GDPR, 'LineWidth', 1.5);
title('Annual Growth Rate of Real GDP');
xlabel('Year');
ylabel('Growth Rate (%)');
grid on;


% MODELO
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% True values 
y    = GDPR(4:end,1);
y1   = GDPR(3:end-1,1);
y2   = GDPR(2:end-2,1);
y3   = GDPR(1:end-3,1);
c    = ones(size(y,1),1);
t    = size(y,1);           % sample length
k    = 3;                   % # RHS variables
x    = [c y1 y2];           % explanatory variable: i.i.d. draws from N(0,1)

tau2 = 1;                   % starting value for Gibbs sampler.
rho  = 0.8;                 % starting value for Gibbs sampler .

% Bayesian Estimation
% Prior for rho
h0 = 0.5;     % mean of prior distribution for rho
M0 = 1;       % variance of prior distribution for rho

% Prior for beta

b0 = [1 0.6 0];                % mean of prior distribution for beta
P0 = [1 0 0; 0 1 0; 0 0 1];       % variance of prior distribution for beta

% Prior for tau2

t0 = 0;               % prior shape parameter for tau
R0 = 0;               % prior scale parameter for tau: zero implies a "flat" prior 
                      % i.e. no prior information<

% Set parameters for Gibbs sampler
nburn = 900000;       % burn-in period: draws are discarded
nkeep = 100000;       % numbers of draws on which inference will be based
n     = nkeep+nburn;  % total number of Gibbs sampling iterations

% Store 
Beta  = zeros(nkeep,k);
Tau2  = zeros(nkeep,1);
Rho   = zeros(nkeep,1);

% Gibb Sampling
for iter = 1:n
    % We clear the error in the AR(2) model and we introduce it into the
    % correlated errors function, we obtain these new RHS variables
    y_new  = y - rho*y1;
    c_new  = c - rho;
    y1_new = y1 - rho*y2;
    y2_new = y2 - rho*y3;
    X_new  = [c_new y1_new y2_new];
    
    % Generate a draw for beta conditional on sig2, the posterior
    P1     = inv(inv(P0) + (inv(tau2))*(X_new'*X_new));
    b1     = (P1)*((inv(P0))*b0' + (inv(tau2))*(X_new'*y_new));
   
% How to take draws in normal distribution with k=3
    beta   = b1 + chol(P1)*randn(k,1);
%         if beta<1 
%            chk=1; 
%         end
%     end
    e_new2 = y - x*beta;
    e_new  = e_new2(2:end,1);
    E_new  = e_new2(1:end-1,1);

% Generate a draw for rho conditional on beta and sig2
    M1  = inv(inv(M0) + (inv(tau2))*(E_new'*E_new)) ;
    h1  = (M1)*((inv(M0))*h0 + (inv(tau2))*(E_new'*e_new));   

% How to take draws in normal distribution    
    rho = h1 + chol(M1)*randn(1,1);
    
% Generate a draw for sig2 conditional on beta
    t1      = t0+t;
    R1      = R0+(y_new-X_new*beta)'*(y_new-X_new*beta);
    shpe    = t1;
    scle    = inv(R1);  
    tau2inv = gamrnd(shpe,scle);
    tau2    = 1/tau2inv;
    
    % Saving data of beta, rho and tau2
    if iter>nburn
        Beta(iter-nburn,:) = beta;
        Rho(iter-nburn,:)  = rho;
        Tau2(iter-nburn,:) = tau2;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% B.2. Histograms for beta, rho and tau^2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(2)
dim1  = 2;
dim2  = 3;
nbins = 500;

%Plot for rho
subplot(dim1,dim2,1)
hist(Rho,nbins);
title('Posterior distribution of \rho')
v=axis;
line(median(Rho)*ones(1,2),[0 max(hist(Rho,nbins))],'Linewidth',1,'Color','r')
line(h0*ones(1,2),[0 max(hist(Rho,nbins))],'Linewidth',1,'Color','y')
legend('estimated','posterior mean','prior')
xlabel('\rho value')
ylabel('Frequency')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','#0072BD','EdgeColor','#0072BD')

% Plot for constant
subplot(dim1,dim2,2)
hist(Beta(:,1),nbins);
title('Posterior distribution of constant')
v=axis;
line(median(Beta(:,1))*ones(1,2),[0 max(hist(Beta(:,1),nbins))],'Linewidth',1,'Color','r')
line(b0(:,1)*ones(1,2),[0 max(hist(Beta(:,1),nbins))],'Linewidth',1,'Color','y')
legend('estimated','posterior mean','prior')
xlabel('b_1 value')
ylabel('Frequency')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','#0072BD','EdgeColor','#0072BD')

%Plot for b_1
subplot(dim1,dim2,3)
hist(Beta(:,2),nbins);
title('Posterior distribution of b_1')
v=axis;
line(median(Beta(:,2))*ones(1,2),[0 max(hist(Beta(:,2),nbins))],'Linewidth',1,'Color','r')
line(b0(:,2)*ones(1,2),[0 max(hist(Beta(:,2),nbins))],'Linewidth',1,'Color','y')
legend('estimated','posterior mean','prior')
xlabel('b_2 value')
ylabel('Frequency')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','#0072BD','EdgeColor','#0072BD')

%Plot for b_2
subplot(dim1,dim2,4)
hist(Beta(:,3),nbins);
title('Posterior distribution of b_2')
v=axis;
line(median(Beta(:,3))*ones(1,2),[0 max(hist(Beta(:,3),nbins))],'Linewidth',1,'Color','r')
line(b0(:,3)*ones(1,2),[0 max(hist(Beta(:,3),nbins))],'Linewidth',1,'Color','y')
legend('estimated','posterior mean','prior')
xlabel('b_3 value')
ylabel('Frequency')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','#0072BD','EdgeColor','#0072BD')

% Plot for tau^2
subplot(dim1,dim2,5)
hist(Tau2,nbins);
title('Posterior distribution of \tau^2')
line(median(Tau2)*ones(1,2),[0 max(hist(Tau2,nbins))],'Linewidth',1,'Color','r')
line(tau2*ones(1,2),[0 max(hist(Tau2,nbins))],'Linewidth',1,'Color','#77AC30')
legend('estimated','posterior mean','true')
xlabel('\tau^2 value')
ylabel('Frequency')
h = findobj(gca,'Type','patch');
set(h,'FaceColor','#0072BD','EdgeColor','#0072BD')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





















