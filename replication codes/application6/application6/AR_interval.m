function [hupper,hlower] = AR_interval(alpha_iz,alpha_1z,sigma_ii,sigma_i1,...
    sigma_11,alpha_size)
% This function calculates Anderson-Rubin confidence interval
chi2 = chi2inv(1-alpha_size,1);
a = alpha_1z^2 - chi2*sigma_11;
b = -2*alpha_iz*alpha_1z + chi2*2*sigma_i1;
c = alpha_iz^2 - chi2*sigma_ii;
if a <= 0
    'Anderson-Rubin confidence interval not defined'
    stop
end
xb = sqrt(b^2 - 4*a*c);
hupper = (-b +xb)/(2*a);
hlower = (-b - xb)/(2*a);
end
