%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GROUP 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
clc;
load("data2.mat");
data1 = [gdp cpi dff m2];
% Seed in order to replicate the same results
seednumber=372398;
rand('seed',seednumber);
randn('seed',seednumber);

% Renaming the Variables
% GDP = y(:,1);
% INF = y(:,2);
% EFR = y(:,3);
% M2  = y(:,4);
% yy  = data;
time = (1985.25:0.25:2007)';

%% 2.6 Suppose you wanted to identify the shocks underlying this model by means of sign restrictions

    % a. Provide a plot for the impact effect of a one-standard deviation shock using the analytical
      %expression for the implicit prior distribution.

      const=ones(size(data1,1),1);   %deterministic term: constant
      ndet=1;                       %number of deterministic variables
      data=[const data1];
     
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


n = 4;
figure;
tiledlayout(n, n);
P = chol(vmat, 'lower');

shocks = ["Demand", "Supply", "Monetary Policy", "Liquidity"];
variables = ["GDP", "Inflation", "Interest Rate", "M2"];

for j = 1:n
    for i = 1:n
        p_ii = P(i, i);

        h_vals = linspace(-p_ii, p_ii, 500);
        q_vals = h_vals / p_ii;
        
        gamma_n2 = gamma(n/2);
        gamma_12 = gamma(1/2);
        gamma_n12 = gamma((n-1)/2);
        density_q = (gamma_n2 ./ (gamma_12 .* gamma_n12)) .* (1 - q_vals.^2).^((n-3)/2);

        density_h = density_q / p_ii;

        nexttile;
        plot(h_vals, density_h, 'b', 'LineWidth', 2);
        
        grid on;
        xlabel(''); 
        ylabel(variables(i)); 
        title([shocks(j) + " shock"]);
    end
end


folder_name = 'graphs 4';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

file_name = fullfile(folder_name, 'impact_effect.png');
saveas(gcf, file_name);


   
    %% b.Verify empirically what the impact effect for each variable looks like. Report plots of the
      % impact effects and provide the numerical values for the cut-off points.

            const=ones(size(data1,1),1);   %deterministic term: constant
      ndet=1;                       %number of deterministic variables
      data=[const data1];

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
    if ~exist('graphs 4', 'dir')
        mkdir('graphs 4');
    end

    saveas(gcf, 'graphs 4/question 6a IRF demand shock.png');

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
    if ~exist('graphs 3', 'dir')
        mkdir('graphs 3');
    end

    saveas(gcf, 'graphs 4/question 6a IRF supply shock .png');
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
    if ~exist('graphs 3', 'dir')
        mkdir('graphs 3');
    end

    saveas(gcf, 'graphs 4/question 6a IRF monetary policy .png');
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
        if ~exist('graphs 3', 'dir')
        mkdir('graphs 3');
    end

    saveas(gcf, 'graphs 4/question 6a IRF money stock .png');
    
    vmat
shocks = ["Demand", "Supply", "Monetary Policy", "Liquidity"];
variables = ["GDP", "Inflation", "FFR", "M2"];
stats_labels = ["Min", "Max", "Median"];

stats = zeros(12, 4);

row_names = strings(12, 1);
row_index = 1;

for i = 1:4  
    for s = 1:3  
        row_names(row_index) = variables(i) + " - " + stats_labels(s);
        for j = 1:4  
            data = squeeze(irf(i, j, 1, :)); 
            if s == 1
                stats(row_index, j) = min(data);
            elseif s == 2
                stats(row_index, j) = max(data);
            elseif s == 3
                stats(row_index, j) = median(data);
            end
        end
        row_index = row_index + 1;
    end
end

table_data = array2table(stats, ...
    'VariableNames', shocks, ...
    'RowNames', row_names);

if ~exist('graphs 4', 'dir')
    mkdir(['graphs 4']);
end

file_name = fullfile('graphs 4', 'IRF_Statistics.xlsx');
writetable(table_data, file_name, 'WriteRowNames', true);
