% Definir parámetros
a = 0.2;
phi = 1;
delta = 0.05;
rho = 0.01;
sigma_ss = 0.2;
b = 0.05;
nu = 0.02;
sigma = 1; % Asegúrate de definir sigma antes de usarlo

% Calcular valores de mu y sig
mu = b * (sigma_ss - sigma);
sig = (1/2) * nu^2 * sigma;

x = linspace(0.001, 0.999, 100) % Cambia este valor según sea necesario
%llamar a la función M
M = buildM(x, mu, sig);

u = sigma^2 * (1 - theta) * theta;

rho*theta=u+M*theta



