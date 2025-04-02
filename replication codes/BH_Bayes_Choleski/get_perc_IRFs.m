IRF_sum = zeros(3,3,hmax,5);

x1=-sort(cumsum(squeeze(IRF(1,1,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(1,1,:,1:3) = temp1;
IRF_sum(1,1,:,4) = temp2;
IRF_sum(1,1,:,5) = temp3;

x1=sort(cumsum(squeeze(IRF(1,2,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(1,2,:,1:3) = temp1;
IRF_sum(1,2,:,4) = temp2;
IRF_sum(1,2,:,5) = temp3;

x1=sort(cumsum(squeeze(IRF(1,3,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(1,3,:,1:3) = temp1;
IRF_sum(1,3,:,4) = temp2;
IRF_sum(1,3,:,5) = temp3;

x1=-sort(cumsum(squeeze(IRF(2,1,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(2,1,:,1:3) = temp1;
IRF_sum(2,1,:,4) = temp2;
IRF_sum(2,1,:,5) = temp3;

x1=sort(cumsum(squeeze(IRF(2,2,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(2,2,:,1:3) = temp1;
IRF_sum(2,2,:,4) = temp2;
IRF_sum(2,2,:,5) = temp3;

x1=sort(cumsum(squeeze(IRF(2,3,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(2,3,:,1:3) = temp1;
IRF_sum(2,3,:,4) = temp2;
IRF_sum(2,3,:,5) = temp3;

x1=-sort(cumsum(squeeze(IRF(3,1,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(3,1,:,1:3) = temp1;
IRF_sum(3,1,:,4) = temp2;
IRF_sum(3,1,:,5) = temp3;

x1=sort(cumsum(squeeze(IRF(3,2,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(3,2,:,1:3) = temp1;
IRF_sum(3,2,:,4) = temp2;
IRF_sum(3,2,:,5) = temp3;

x1=sort(cumsum(squeeze(IRF(3,3,:,:)),1),2);
temp1=[(median(x1,2)) x1(:,index68(1)) x1(:,index68(2))];
temp2 = x1(:,index95(1));
temp3 = x1(:,index95(2));
IRF_sum(3,3,:,1:3) = temp1;
IRF_sum(3,3,:,4) = temp2;
IRF_sum(3,3,:,5) = temp3;

clear IRF

IRF = IRF_sum;