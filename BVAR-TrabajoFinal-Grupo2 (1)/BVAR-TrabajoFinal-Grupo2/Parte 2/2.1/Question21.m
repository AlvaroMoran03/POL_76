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
    V1 = getFredData('GDPC1', '1985-04-01', '2006-12-31','pch','q');
    V2 = getFredData('GDPDEF', '1985-04-01', '2006-12-31','pch','q');
    V3 = getFredData('M2REAL', '1985-04-01', '2006-12-31','pch','q');
    V4 = getFredData('FEDFUNDS', '1985-04-01', '2006-12-31','lin','q','avg');

    % Renaming the Variables

    GDP = V1.Data(:,2); DEF = V2.Data(:,2); M2 = V3.Data(:,2); FED = V4.Data(:,2);
    
    clear V1 V2 V3 V4
    
    data=[GDP DEF FED M2];
    
    % Plot the transformed data
    figure(1)
    plot_graphs(data(:,1),data(:,2),data(:,2),data(:,4));

    % b. Write Matlab code to estimate a reduced-form VAR(4) model with a constant term using OLS. Report the
    % estimates of the reduced-form covariance matrix and the (𝑘𝑘×𝑛𝑛) matrix of reduced-form coefficients where
    % nn=4 the number of endogenous variables.

    % Preparing the data
    yy  = data(5:end, :);
    T   = size(yy,1);
    xx  = [ones(T,1) data(4:end - 1, :) data(3:end - 2, :) data(2:end - 3, :) data(1:end - 4, :)];
    k   = size(xx,2);

    % Estimating the VAR(4) by OLS
    b_ols = inv(xx'*xx)*xx'*yy;
    resid = yy - xx*b_ols;
    omega_ols = (resid'*resid)/(T-k);
    