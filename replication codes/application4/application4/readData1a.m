nlags = 8;  % this is number of lags in VAR
kcum = 1;   % This will accumulate IRFs and var decomps
load labor_data.csv;   % 1947:Q1 to 2014:Q2
% col 1 = date 
% col 2 = real labor compensation 
% col 3 = total employment
wage = 100*(log(labor_data(2:end,2)) - log(labor_data(1:end-1,2)));
employment = 100*(log(labor_data(2:end,3)) - log(labor_data(1:end-1,3)));

varnames = {' wage'; ' employment'};
shocknames = {' demand'; ' supply'};
yall = [wage employment];
tstart = 92;    % start estimation with 1970:Q1
tend = size(yall,1);    % end estimation with 2014:Q2
YY = yall(tstart:tend,:);
date_start = 1970;
date_end = 2014.25;
time = (date_start:0.25:date_end)';
