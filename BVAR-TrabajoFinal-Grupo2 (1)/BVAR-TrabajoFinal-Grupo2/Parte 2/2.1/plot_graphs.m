function [] = plot_graphs(y1,y2,y3, y4)

%date
t1 = datetime(1959,1,1);
t2 = datetime(2006,4,1);
date = t1 + calquarters(1:47.8*4);

ia=2;
ib=2;
%tt=1:16;

subplot(ia,ib,1)
plot(date,y1,'.-')
%legend('contemporaneous shock','8 - quarter in advance signal')
axis tight
ylabel('Quarterly variation');
set (gca,'XTicklabel',{''})
title('US real GDP');

subplot(ia,ib,2)
plot(date,y2,'.-');
axis tight
ylabel('Quarterly variation');
set (gca,'XTicklabel',{''})
title('Inflation');

subplot(ia,ib,3)
plot(date,y3,'.-');
axis tight
ylabel('Percent');
set (gca,'XTicklabel',{''})
title('Federal Funds Rate');

subplot(ia,ib,4)
plot(date,y4,'.-');
axis tight
ylabel('Quarterly variation');
set (gca,'XTicklabel',{''})
title('M2 money stock');

end

