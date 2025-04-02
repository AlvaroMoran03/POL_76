%% Paso 1: Definir parámetros
rng(1234)
a = 0.2;
phi = 1;
delta = 0.05;
rho = 0.01;
sigma_ss = 0.2;
b = 0.05;
nu = 0.02;

%% Paso 2: Crear grid para sigma (ahora `x` representa `sigma`)
N = 1000;
x = linspace(0.0001, 0.9999, N)'; % Grid entre 0.001 y 0.999

%% Paso 3: Construir la matriz M usando buildM
mu_o = b * (sigma_ss - x); % Drift term
sig_o = nu * (x.^(0.5)); % Diffusion term

M = buildM(x, mu_o, sig_o); % Construcción de la matriz M

%% Paso 4: Definir u(θ) basado en la ecuación encontrada en (b)
u_theta = @(theta, sigma) sigma.^2 .* ((1 - theta).^2) .* theta; % Función de utilidad

%% Paso 5: Configurar la iteración de la función de valor con dt decreciente
dt_0 = 0.1;   % Paso de tiempo inicial
dt = dt_0;    % Inicialización de dt
tol = 1e-6;   % Tolerancia de convergencia
max_iter = 1000; % Máximo de iteraciones
decay_factor = 0.99; % Factor de reducción de dt

%% Definir theta inicial aleatorio
theta = rand(N, 1);

%% Iteración de la ecuación de Bellman con dt decreciente
for iter = 1:max_iter
    % Disminuir dt progresivamente
    dt = dt * decay_factor;  % Reducción en cada iteración

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

%% Mostrar resultado final
disp('Theta final después de iteraciones:');
disp(theta);

%% Graficar theta(σ̃) si deseas visualizar el resultado
figure;
plot(x, theta, 'LineWidth', 2);
xlabel('\sigmã');
ylabel('\vartheta(\sigmã)');
title('Solución de \vartheta(\sigmã) con Value Function Iteration y dt fijo');
grid on;   

