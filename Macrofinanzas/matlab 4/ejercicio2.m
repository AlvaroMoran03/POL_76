clear; clc;

%% 🔹 Definir valores de los parámetros
rho_e = 0.06;  
rho_h = 0.05; 
a_e   = 0.11; 
a_h   = 0.03;  
delta = 0.05;  
phi   = 10;    
alpha = 0.5;   
sigma = 0.1;  

%% 🔹 Definir grilla uniforme para eta
N   = 1000;
eta_vals = linspace(0.0001, 0.9999, N)';
q_vals = zeros(size(eta_vals));

% Condición inicial de q
q_vals(1) = (1 + a_h * phi) / (1 + rho_h * phi);

% Definir la función iota
iota = @(q) (1/phi) * q - (1/phi);


kappa_vals = zeros(size(eta_vals));

for i = 2:length(eta_vals)
    eta = eta_vals(i);
    eta_prev = eta_vals(i-1);
    q_prev = q_vals(i-1);
    
    % Aproximación de la derivada dq/deta
    dq_deta = @(q) (q - q_prev) / (eta - eta_prev);
    
    % Definir la función kappa
    kappa = @(q) (iota(q) + q * (eta * rho_e + (1 - eta) * rho_h) - a_h) / (a_e - a_h);
    
    % Definir la ecuación a resolver
    func = @(q) ((a_e - a_h) / q - ((kappa(q) - eta) / (eta * (1 - eta))) * ((sigma / (1 - ((kappa(q) / eta) - 1) * (eta * dq_deta(q) / q)))^2));
    
    % Resolver para q usando fsolve
    options = optimoptions('fsolve', 'Display', 'off');
    q_vals(i) = fsolve(func, q_prev, options);
    
    % Guardar el valor de kappa para q encontrado
    kappa_vals(i) = kappa(q_vals(i));
end

% Graficar la solución
figure;
plot(eta_vals, q_vals, 'b', 'LineWidth', 2);
xlabel('\eta');
ylabel('q');
title('Solución de q para cada \eta');
grid on;



