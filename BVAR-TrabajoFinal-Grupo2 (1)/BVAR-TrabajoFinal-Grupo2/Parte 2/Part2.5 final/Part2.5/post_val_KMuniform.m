function objective=post_val_KMuniform(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M4,M_star1,M_star2,M_star3,M_star4,S)

 a_alpha=param(1);
 b_alpha=param(2);

 T=param(3);
 omega_hatT=reshape(param(4:19),4,4);
 kappastar=param(20);
 kappa=param(21);

 prior_a = 1/(b_alpha-a_alpha);
 prior_j = 1/(b_alpha-a_alpha);
 prior_k = 1/(b_alpha-a_alpha);
 prior_l = 1/(b_alpha-a_alpha);
 prior_e = 1/(b_alpha-a_alpha);
 prior_b = 1/(b_alpha-a_alpha);
 prior_c = 1/(b_alpha-a_alpha);
 prior_d = 1/(b_alpha-a_alpha);
 prior_f = 1/(b_alpha-a_alpha);
 prior_g = 1/(b_alpha-a_alpha);
 prior_h = 1/(b_alpha-a_alpha);
 prior_i = 1/(b_alpha-a_alpha);
%A=[1 0 -theta_hat(1,1); 0 1 -theta_hat(2,1); -theta_hat(3:4,1)' 1];
A=[1 -theta_hat(1:3,1)'; -theta_hat(4,1) 1 -theta_hat(5:6,1)'; -theta_hat(7:8,1)' 1 -theta_hat(9,1); -theta_hat(10:12,1)' 1 ];

omega=A*S*A';
taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);

log_priors=log(prior_a)+log(prior_b)+log(prior_c)+log(prior_d)+log(prior_e)+log(prior_f)+log(prior_g)+log(prior_h)+log(prior_i)+log(prior_j)+log(prior_k)+log(prior_l);

up=log_priors+T/2*log(det(A*omega_hatT*A'));
down=kappastar*log((2/T)*taustar1) ...
    +kappastar*log((2/T)*taustar2) ...
    +kappastar*log((2/T)*taustar3) ...
    +kappastar*log((2/T)*taustar4);
objective=-(up-down);

end









