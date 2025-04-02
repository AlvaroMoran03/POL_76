clear; clc;
%% Guardar graficos
save_folder = fullfile(pwd, 'Graficos 2'); 

% Crear la carpeta si no existe
if ~exist(save_folder, 'dir')
    mkdir(save_folder);
end
%% a Definir valores de los parámetros
rho_e = 0.06;  
rho_h = 0.05; 
a_e   = 0.11; 
a_h   = 0.03;  
delta = 0.05;  
phi   = 10;    
alpha = 0.5;   
sigma = 0.1;  

%% 2b Definir grilla uniforme para eta
N   = 1000;
eta_vals = linspace(0.0001, 0.9999, N)';
%% 2 c 
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
    
    % Resolver para q 
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
    
    % Resolver q_vals(i) 
    q_vals(i) = (a_e + (1/phi)) / ((eta * rho_e + (1 - eta) * rho_h) + (1/phi));
    
    % Llenar kappa_vals con 1
    kappa_vals(i) = 1;
end
%% 2d
iota_vals = zeros(size(q_vals));

% Calcular iota(q) 
for i = 2:length(q_vals)
    iota_vals(i) = (1/phi) * q_vals(i) - (1/phi);
end

sigma_q_vals = zeros(size(q_vals)-1);

% Iterar para hallar sigma_q
for i = 2:length(q_vals)
    eta = eta_vals(i);
    q = q_vals(i);
    dq_deta = (q_vals(i) - q_vals(i-1)) / (eta_vals(i) - eta_vals(i-1)); % Aproximación de dq/deta
    
    % Calcular sigma^q
    sigma_q_vals(i) = (sigma / (1 - ((kappa_vals(i) / eta) - 1) * (eta * dq_deta / q))) - sigma;
end


% Graficar la solución de q
figure;
plot(eta_vals, q_vals, 'r', 'LineWidth', 2); % Rojo
xlabel('\eta');
ylabel('q');
grid on;
saveas(gcf, fullfile(save_folder, 'q.png'));

% Graficar la solución de kappa^e
figure;
plot(eta_vals, kappa_vals, 'g', 'LineWidth', 2);
xlabel('\eta');
ylabel('\kappa^e');
grid on;
saveas(gcf, fullfile(save_folder, 'Kappa_e.png'));

% Graficar la solución de iota
figure;
plot(eta_vals, iota_vals, 'm.', 'LineWidth', 2);
xlabel('\eta');
ylabel('\iota');
grid on;
saveas(gcf, fullfile(save_folder, 'iota.png'));

% Graficar la solución de sigma^q, 
figure;
plot(eta_vals(2:end), sigma_q_vals(2:end), 'b', 'LineWidth', 2); 
xlabel('\eta');
ylabel('\sigma^q');
grid on;
saveas(gcf, fullfile(save_folder, 'sigma_q.png'));

%% f

eta_mu_eta_vals = zeros(size(eta_vals));
eta_sigma_eta_vals = zeros(size(eta_vals));

% Calcular eta * mu_eta y eta * sigma_eta
for i = 2:(length(eta_vals))
    eta = eta_vals(i);
    sigma_q = sigma_q_vals(i);
    kappa_e = kappa_vals(i);  
    
    % Cálculo de eta * mu_eta
    eta_mu_eta_vals(i) = eta * (1 - eta) * (rho_h - rho_e) + ...
        (kappa_e - 2 * eta * kappa_e + eta^2) * ...
        ((kappa_e - eta) / (eta * (1 - eta))) * (sigma + sigma_q)^2;
    
    % Cálculo de eta * sigma_eta
    eta_sigma_eta_vals(i) = (kappa_vals(i) - eta) * (sigma + sigma_q);
end

% Graficar eta * mu_eta
figure;
plot(eta_vals(2:end), eta_mu_eta_vals(2:end), 'r', 'LineWidth', 2);
xlabel('\eta');
ylabel('\eta \mu_\eta');
grid on;
saveas(gcf, fullfile(save_folder, 'eta_mu_eta_plot.png'));

% Graficar eta * sigma_eta
figure;
plot(eta_vals(2:end), eta_sigma_eta_vals(2:end), 'b', 'LineWidth', 2);
xlabel('\eta');
ylabel('\eta \sigma_\eta');
grid on;
saveas(gcf, fullfile(save_folder, 'eta_sigma_eta_vals.png'));

%% g


% Inicializar el vector para almacenar los valores de r
r_vals = zeros(size(eta_vals)-2);

% Aproximar derivadas y calcular r
for i = 2:length(eta_vals)-1
    % Extraer valores para el índice i
    eta = eta_vals(i);
    q = q_vals(i);
    kappa = kappa_vals(i);
    iota_q = iota_vals(i);
    sigma_q = sigma_q_vals(i);
    eta_mu_eta = eta_mu_eta_vals(i);
    eta_sigma_eta = eta_sigma_eta_vals(i);
    
    % Valores en los puntos vecinos
    eta_prev = eta_vals(i-1);
    eta_next = eta_vals(i+1);
    q_prev = q_vals(i-1);
    q_next = q_vals(i+1);
    
    % Aproximación de la primera derivada q_eta usando diferencia centrada
    q_eta = (q - q_prev) / (eta - eta_prev);
    
    % Aproximación de la segunda derivada q_etaeta usando diferencias finitas
    q_etaeta = (q_next - 2*q + q_prev) / ((eta_next- eta_prev)^2);
    
    % Definir Phi(iota) según la ecuación dada
    Phi_iota = (1 / phi) * log(1 + phi * iota_q);
    
    % Calcular r según la ecuación dada
    r_vals(i) = (q_eta/q) * (eta_mu_eta) + ...
                (1 / (2*q)) * q_etaeta * (eta_sigma_eta)^2 + ...
                Phi_iota + sigma * sigma_q + ...
                (kappa * a_e + (1 - kappa) * a_h) / q - ...
                (((kappa^2 / eta) + ((1 - kappa)^2 / (1 - eta))) * (sigma + sigma_q)^2);
end

figure;
plot(eta_vals(2:end), r_vals, 'r', 'LineWidth', 2); % Azul con línea punteada
xlabel('\eta');
ylabel('r');
grid on;
saveas(gcf, fullfile(save_folder, 'r.png'));
%Sacar outliers
Q1 = quantile(r_vals(2:end), 0.25); 
Q3 = quantile(r_vals(2:end), 0.75); 
IQR = Q3 - Q1; 

% Definir límites para considerar valores no outliers
lower_bound = Q1 - 1.5 * IQR;
upper_bound = Q3 + 1.5 * IQR;

% Filtrar valores dentro del rango permitido
valid_indices = (r_vals(2:end) >= lower_bound) & (r_vals(2:end) <= upper_bound);
filtered_eta_vals = eta_vals(2:end);
filtered_r_vals = r_vals(2:end);

% Aplicar filtro de outliers
filtered_eta_vals = filtered_eta_vals(valid_indices);
filtered_r_vals = filtered_r_vals(valid_indices);

% Graficar sin outliers
figure;
plot(filtered_eta_vals, filtered_r_vals, 'r', 'LineWidth', 2);
xlabel('\eta');
ylabel('r');
grid on;

% Guardar el gráfico
saveas(gcf, fullfile(save_folder, 'r_sin_outliers.png'));