function objective=post_val_KMuniform(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M4,M_star1,M_star2,M_star3,M_star4,S)

c_general=param(1);
sigma_general=param(2);
nu_general=param(3);
c_general=param(4);
sigma_general=param(5);
nu_general=param(6); 
c_general=param(7);
sigma_general=param(8);
nu_general=param(9);
c_general=param(10);
sigma_general=param(11);
nu_general=param(12);
c_general=param(13);
sigma_general=param(14);
nu_general=param(15);
c_general=param(16);
sigma_general=param(17);
nu_general=param(18);
T=param(19);
omega_hatT=reshape(param(20:35),4,4);
kappastar=param(36);
kappa=param(37);

%prior_a = 1/(b_alpha-a_alpha);
prior_r_m_1 = student_pos_prior(theta_hat(1,1),c_general,sigma_general,nu_general);
prior_r_pi_1 = student_pos_prior(theta_hat(2,1),c_general,sigma_general,nu_general);
prior_r_y_1 = student_pos_prior(theta_hat(3,1),c_general,sigma_general,nu_general);
prior_m_r_1 = student_neg_prior(theta_hat(4,1),c_general,sigma_general,nu_general);
prior_m_pi_1 = student_pos_prior(theta_hat(5,1),c_general,sigma_general,nu_general);
prior_m_y_1 = student_pos_prior(theta_hat(6,1),c_general,sigma_general,nu_general);
prior_pi_r_1 = student_neg_prior(theta_hat(7,1),c_general,sigma_general,nu_general);
prior_pi_m_1 = student_neg_prior(theta_hat(8,1),c_general,sigma_general,nu_general);
prior_pi_y_1 = student_neg_prior(theta_hat(9,1),c_general,sigma_general,nu_general);
prior_y_r_1 = student_neg_prior(theta_hat(10,1),c_general,sigma_general,nu_general);
prior_y_m_1 = student_neg_prior(theta_hat(11,1),c_general,sigma_general,nu_general);
prior_y_pi_1 = student_pos_prior(theta_hat(12,1),c_general,sigma_general,nu_general);

    
A_1=[1 -theta_hat(1:3,1)'; -theta_hat(4,1) 1 -theta_hat(5:6,1)'; -theta_hat(7:8,1)' 1 -theta_hat(9,1); -theta_hat(10:12,1)' 1 ];

omega=A_1*S*A_1';
taustar1=gettau(kappa,omega(1,1),A_1(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A_1(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A_1(3,:),Syy,Sxy,eta,M3,M_star3);
taustar4=gettau(kappa,omega(4,4),A_1(4,:),Syy,Sxy,eta,M4,M_star4);

log_priors= log(prior_r_m_1)+log(prior_r_pi_1)+log(prior_r_y_1)+log(prior_m_pi_1)+log(prior_m_y_1)+log(prior_pi_y_1)+log(prior_m_r_1)+log(prior_pi_r_1)+log(prior_pi_m_1)+log(prior_y_r_1)+log(prior_y_m_1)+log(prior_y_pi_1);
up=log_priors+T/2*log(det(A_1*omega_hatT*A_1'));
down=kappastar*log((2/T)*taustar1) ...
    +kappastar*log((2/T)*taustar2) ...
    +kappastar*log((2/T)*taustar3) ...
    +kappastar*log((2/T)*taustar4);
objective=-(up-down);











