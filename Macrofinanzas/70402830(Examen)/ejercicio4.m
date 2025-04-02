%% a Definir parámetros
rng(1234)
a = 0.2;
phi = 1;
delta = 0.05;
rho = 0.01;
sigma_ss = 0.2;
b = 0.05;
nu = 0.02;

%% c-i 
% Crear grid para sigma (ahora `x` representa `sigma`)
N = 1000;
x = linspace(0.0001, 0.9999, N)'; % Grid entre 0.001 y 0.999

 %Construir la matriz M usando buildM
mu_o = b * (sigma_ss - x); % Drift term
sig_o = nu * (x.^(0.5)); % Diffusion term

M = buildM(x, mu_o, sig_o); % Construcción de la matriz M

%% c-iii
u_theta = @(theta, sigma) sigma.^2 .* ((1 - theta).^2) .* theta; % Función de utilidad

%% c iii-iv
dt = 0.1;   
tol = 1e-6;   
max_iter = 1000; 
%definir un theta aleatorio
theta = rand(N, 1);

% Iteración de la ecuación de Bellman con dt decreciente
for iter = 1:max_iter

    % Matriz del método implícito con dt variable
    I = eye(N);
    A = (1 + rho * dt) * I - dt * M; % Matriz a invertir
    
    % Verificar condición de A antes de resolver
    if rcond(A) < 1e-12
        warning('La matriz A puede estar mal condicionada o ser singular en iteración %d.', iter);
        break;
    end

    % Evaluar u(theta, sigma) en cada iteración
    u_eval = u_theta(theta, x);
    
    % Asegurar que u_eval es un vector columna
    u_eval = u_eval(:);
    
    % Resolver para theta en esta iteración
    theta_new = A \ (dt * u_eval + theta);
    
    % Criterio de convergencia
    diff = norm(theta_new - theta, Inf);
    if diff < tol
        fprintf('Convergencia alcanzada en %d iteraciones con dt = %.6f.\n', iter, dt);
        break;
    end
    
    % Actualizar theta sin mezcla
    theta = theta_new;
end



% Crear carpeta "Graficos 4" si no existe
folder_name = 'Graficos 4';
if ~exist(folder_name, 'dir')
    mkdir(folder_name);
end

%% d Graficar 
% theta(σ̃) 
figure;
plot(x(1:970), theta(1:970), 'LineWidth', 2);
xlabel('\sigmã');
ylabel('\vartheta(\sigmã)');
grid on;   
saveas(gcf, fullfile(folder_name, 'theta_plot.png'));

% Calcular q^k y q_t^B para cada valor de vartheta
q_k_vals = (1 - theta) .* (a * phi + 1) ./ (phi * rho + 1 - theta);
q_B_vals = (theta .* (a * phi + 1)) ./ (phi * rho + 1 - theta);

% Graficar q^k 
figure;
plot(x (1:970), q_k_vals(1:970), 'b', 'LineWidth', 2);
xlabel('\sigmã');
ylabel('q^k');
grid on;
saveas(gcf, fullfile(folder_name, 'q_k_plot.png')); 

% Graficar q_t^B 
figure;
plot(x(1:970), q_B_vals(1:970), 'r', 'LineWidth', 2);
xlabel('\sigmã');
ylabel('q^B');
grid on;
saveas(gcf, fullfile(folder_name, 'q_t_B_plot.png')); 
% Graficar z_tilde

zeta_tilde = (1 - theta).* x;

figure;
plot(x(1:970), zeta_tilde(1:970), 'r', 'LineWidth', 2);
xlabel('\sigmã');
ylabel('$\tilde{\zeta}_t$', 'Interpreter', 'latex');
grid on;
saveas(gcf, fullfile(folder_name, 'zeta.png')); 
% Graficar r_f 

rf = rho + (1 / phi) * log(1 + phi * ((1 / phi) * q_k_vals - (1 / phi))) -delta - ((1 - theta) .*x).^2;

figure;
plot(x (1:970), rf (1:970), 'b', 'LineWidth', 2);
xlabel('\sigmã');
ylabel('r^f');
grid on;
saveas(gcf, fullfile(folder_name, 'rf.png')); 
