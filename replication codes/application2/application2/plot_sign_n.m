clear;
clc;

grid_size = 0.01;
xgrid = (-1+grid_size):grid_size:(1 - grid_size);
n=2
%Implied Distribution for q
y2 = (1 - (xgrid.^2)).^((n-3)/2);
a=gamma(n/2)/(gamma(1/2)*gamma((n-1)/2));
y2 = y2*a;

subplot(2,1,1)
plot(xgrid,y2,'b-','LineWidth',1.5)
set(gca,'Xtick',0)
text(-1 ,-0.25,'-\sigma_{ii}^{0.5}','hor','center')
text(1 ,-0.25,'\sigma_{ii}^{0.5}','hor','center')
hold on

grid_size = 0.01;
xgrid = (-1+grid_size):grid_size:(1 - grid_size);
n=6
y2 = (1 - (xgrid.^2)).^((n-3)/2);
a=gamma(n/2)/(gamma(1/2)*gamma((n-1)/2));
y2 = y2*a;

plot(xgrid,y2,'k--','LineWidth',1.5)
hold on

grid_size = 0.01;
xgrid = (-1+grid_size):grid_size:(1 - grid_size);
n=12
y2 = (1 - (xgrid.^2)).^((n-3)/2);
a=gamma(n/2)/(gamma(1/2)*gamma((n-1)/2));
y2 = y2*a;
plot(xgrid,y2,'r:','LineWidth',1.5)

legend('#variables=2','#variables=6','#variables=12','Location','BestOutside')
title('Response of variable i to any 1SD structural shock')