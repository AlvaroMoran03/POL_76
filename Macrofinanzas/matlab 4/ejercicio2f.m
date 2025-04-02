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
eta_vals(1)=0
% Definir la función iota
iota = @(q) (1/phi) * q - (1/phi);

% Inicializar el vector para almacenar los valores de kappa
kappa_vals = zeros(size(eta_vals));

for i = 2:min(400, length(eta_vals))
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

    % Si kappa es mayor o igual a 1, imponer 1 y detener la simulación
    if kappa_vals(i) >= 1
        kappa_vals(i) = 1;
        break;
    end
end
for i = 401:length(eta_vals)
    eta = eta_vals(i);
    
    % Resolver q_vals(i) despejando de la ecuación dada
    q_vals(i) = (a_e + (1/phi)) / ((eta * rho_e + (1 - eta) * rho_h) + (1/phi));
    
    % Llenar kappa_vals con 1
    kappa_vals(i) = 1;
end

iota_vals = zeros(size(q_vals));

% Calcular iota(q) para cada q_vals
for i = 1:length(q_vals)
    iota_vals(i) = (1/phi) * q_vals(i) - (1/phi);
end

sigma_q_vals = zeros(size(q_vals)-2);

% Iterar desde el segundo hasta el ultimo
for i = 2:length(q_vals)
    eta = eta_vals(i);
    q = q_vals(i);
    dq_deta = (q_vals(i) - q_vals(i-1)) / (eta_vals(i) - eta_vals(i-1)); % Aproximación de dq/deta
    
    % Calcular sigma^q
    sigma_q_vals(i-1) = (sigma / (1 - ((kappa_vals(i) / eta) - 1) * (eta * dq_deta / q))) - sigma;
end


% Graficar la solución de q
figure;
plot(eta_vals, q_vals, 'r', 'LineWidth', 2); % Rojo
xlabel('\eta');
ylabel('q');
title('Solución de q para cada \eta');
grid on;

% Graficar la solución de kappa^e
figure;
plot(eta_vals, kappa_vals, 'g', 'LineWidth', 2); % Verde con línea discontinua
xlabel('\eta');
ylabel('\kappa^e');
title('Solución de \kappa^e para cada \eta');
grid on;

% Graficar la solución de iota
figure;
plot(eta_vals, iota_vals, 'm-.', 'LineWidth', 2); % Magenta con línea de puntos y rayas
xlabel('\eta');
ylabel('\iota');
title('Solución de \iota para cada \eta');
grid on;

% Graficar la solución de sigma^q, omitiendo el primer valor de eta
figure;
plot(eta_vals(2:end), sigma_q_vals, 'b', 'LineWidth', 2); % Azul con línea punteada
xlabel('\eta');
ylabel('\sigma^q');
title('Solución de \sigma^q para cada \eta');
grid on;

%% f

% Inicializar los vectores para almacenar los valores calculados
eta_mu_eta_vals = zeros(size(eta_vals)-1);
eta_sigma_eta_vals = zeros(size(eta_vals)-1);

% Calcular eta * mu_eta y eta * sigma_eta
for i = 2:length(eta_vals)  % Omitimos el primer y último valor
    eta = eta_vals(i);
    sigma_q = sigma_q_vals(i-1);
    kappa_e = kappa_vals(i);  % Asumimos kappa^e = kappa_vals(i)
    
    % Cálculo de eta * mu_eta
    eta_mu_eta_vals(i) = eta * (1 - eta) * (rho_h - rho_e) + ...
        (kappa_e - 2 * eta * kappa_e + eta^2) * ...
        ((kappa_e - eta) / (eta * (1 - eta))) * (sigma + sigma_q)^2;
    
    % Cálculo de eta * sigma_eta
    eta_sigma_eta_vals(i) = (kappa_vals(i) - eta) * (sigma + sigma_q);
end

% Graficar eta * mu_eta
figure;
plot(eta_vals(2:end-1), eta_mu_eta_vals(2:end-1), 'r', 'LineWidth', 2);
xlabel('\eta');
ylabel('\eta \mu_\eta');
title('Solución de \eta \mu_\eta para cada \eta');
grid on;

% Graficar eta * sigma_eta
figure;
plot(eta_vals(2:end-1), eta_sigma_eta_vals(2:end-1), 'b', 'LineWidth', 2);
xlabel('\eta');
ylabel('\eta \sigma_\eta');
title('Solución de \eta \sigma_\eta para cada \eta');
grid on;

