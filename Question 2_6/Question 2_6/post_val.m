%%%%%%%%%%%%%%%%%%%%%%%%%%%% Final Exam %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Course: A Bayesian Approach to Identification of Structural VAR Models
% Teacher: Christiane Baumeister
% Group: 6
% Members: Andreu Marquez, Ricardo Suarez, Alonso Palacios, Mayra Gonzales, Christian

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part II %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function objective=post_val(theta_hat,param,Syy,Sxy,eta,M1,M2,M3,M4,M_star1,M_star2,M_star3,M_star4,S)

c_alpha_piy = param(1);
sigma_alpha_piy = param(2);   
nu_alpha_piy = param(3);

c_alpha_my = param(4);
sigma_alpha_my = param(5);   
nu_alpha_my = param(6);

c_alpha_mpi = param(7);
sigma_alpha_mpi = param(8);   
nu_alpha_mpi = param(9);

c_alpha_ry = param(10);
sigma_alpha_ry = param(11);   
nu_alpha_ry = param(12);

c_alpha_rpi = param(13);
sigma_alpha_rpi = param(14);   
nu_alpha_rpi = param(15);

c_alpha_rm = param(16);
sigma_alpha_rm = param(17);   
nu_alpha_rm = param(18);

T=param(19);
omega_hatT=reshape(param(20:35),4,4);
kappastar=param(36);
kappa=param(37);

prior_a = student_prior(theta_hat(1,1),c_alpha_piy,sigma_alpha_piy,nu_alpha_piy);
prior_b = student_prior(theta_hat(2,1),c_alpha_my,sigma_alpha_my,nu_alpha_my);
prior_c = student_prior(theta_hat(3,1),c_alpha_mpi,sigma_alpha_mpi,nu_alpha_mpi);
prior_d = student_prior(theta_hat(4,1),c_alpha_ry,sigma_alpha_ry,nu_alpha_ry);
prior_e = student_prior(theta_hat(5,1),c_alpha_rpi,sigma_alpha_rpi,nu_alpha_rpi);
prior_f = student_prior(theta_hat(6,1),c_alpha_rm,sigma_alpha_rm,nu_alpha_rm);

A=setA(theta_hat); % specified in setA

omega=A*S*A';
taustar1=gettau(kappa,omega(1,1),A(1,:),Syy,Sxy,eta,M1,M_star1);
taustar2=gettau(kappa,omega(2,2),A(2,:),Syy,Sxy,eta,M2,M_star2);
taustar3=gettau(kappa,omega(3,3),A(3,:),Syy,Sxy,eta,M3,M_star3);
taustar4=gettau(kappa,omega(4,4),A(4,:),Syy,Sxy,eta,M4,M_star4);

log_priors = log(prior_a) + log(prior_b) + log(prior_c) + log(prior_d) + log(prior_e) + log(prior_f);

up=log_priors+T/2*log(det(A*omega_hatT*A'));
down=kappastar*log((2/T)*taustar1) ...
    +kappastar*log((2/T)*taustar2) ...
    +kappastar*log((2/T)*taustar3) ...
    +kappastar*log((2/T)*taustar4);
objective=-(up-down);









