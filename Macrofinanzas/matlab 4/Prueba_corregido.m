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
eta = linspace(0.0001, 0.9999, N)';

%% 🔹 Condiciones de frontera para q
q_0 = (1 + a_h * phi) / (1 + rho_h * phi);  
q_1 = (1 + a_e * phi) / (1 + rho_e * phi); 

%% 🔹 Inicializar variables
q = zeros(N, 1);
q(1) = q_0;  
q_eta = zeros(N, 1);
kappa_e_vals = zeros(N, 1); 

%% 🔹 Primera iteración (de 2 a 400)
for i = 2:400
    % Aproximación de q_eta
    q_eta_approx = @(q_val) (q_val - q(i-1)) / (eta(i) - eta(i-1));

    % Definir la función para calcular kappa_e
    kappa_e_fun = @(q_val) min(((q_val - 1) / phi + q_val * (eta(i) * rho_e + (1 - eta(i)) * rho_h) - a_h) / (a_e - a_h), 1);

    % Definir la ecuación a resolver
    fun = @(q_i) ((a_e - a_h) / q_i) - ...
        ((kappa_e_fun(q_i) - eta(i)) / (eta(i) * (1 - eta(i)))) * ...
        ((sigma / (1 - ((kappa_e_fun(q_i) / eta(i) - 1) * (eta(i) * q_eta_approx(q_i) / q_i))))^2);

    % Resolver ecuación mediante fsolve
    options = optimset('Display', 'off');
    q_sol = fsolve(fun, q(i-1), options);
    
    % Actualizar valores
    q(i) = q_sol;
    q_eta(i) = (q(i) - q(i-1)) / (eta(i) - eta(i-1));

    % Calcular kappa_e
    kappa_e_vals(i) = (q(i) - 1) / phi + q(i) * (eta(i) * rho_e + (1 - eta(i)) * rho_h) - a_h;
    kappa_e_vals(i) = kappa_e_vals(i) / (a_e - a_h);
    
    if kappa_e_vals(i) >= 1
        kappa_e_vals(i:end) = 1; 
        q(i:end) = q(i); 
        q_eta(i:end) = 0; 
        break;
    end
end

%% 🔹 Segunda iteración (de 401 a N) resolviendo q explícitamente
for i = 401:N
    q(i) = (phi * a_e+1) ./ (1 + phi * (eta(i) * rho_e + (1 - eta(i)) * rho_h));  % Solución explícita
    q_eta(i) = (q(i) - q(i-1)) / (eta(i) - eta(i-1));
    kappa_e_vals(i) = 1; % Fijar kappa_e en 1 en esta segunda iteración
end

%% 🔹 Calcular variables adicionales
i_val = (q - 1) / phi;  % Cálculo de i
sigma_q = zeros(N, 1); % Inicializar sigma_q

for i = 1:N
    if eta(i) > 0
        sigma_q(i) = (sigma/ (1-(kappa_e_vals(i) / eta(i) - 1) * ((eta(i) * q_eta(i)) / q(i)))) - sigma;
    else
        sigma_q(i) = NaN;  % Evitar errores en la primera iteración
    end
end

%% 🔹 Graficar los valores calculados
figure;
subplot(3,1,1)
plot(eta, q, 'b', 'LineWidth', 1.5);
xlabel('\eta');
ylabel('q(\eta)');
title('Evolución de q(\eta)');
grid on;

subplot(3,1,2)
plot(eta, i_val, 'r', 'LineWidth', 1.5);
xlabel('\eta');
ylabel('i(\eta)');
title('Evolución de i(\eta)');
grid on;

subplot(3,1,3)
plot(eta, kappa_e_vals, 'g', 'LineWidth', 1.5);
xlabel('\eta');
ylabel('\kappa_e(\eta)');
title('Evolución de \kappa_e(\eta)');
grid on;

figure;
plot(eta, sigma_q, 'm', 'LineWidth', 1.5);
xlabel('\eta');
ylabel('\sigma_q(\eta)');
title('Evolución de \sigma_q(\eta)');
grid on;

