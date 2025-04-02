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
mu = b * (sigma_ss - x); % Drift term
sig = nu * sqrt(x); % Diffusion term

M = buildM(x, mu, sig); % Construcción de la matriz M

%% Paso 4: Definir u(θ) basado en la ecuación encontrada en (b)
u_theta = @(theta, sigma) sigma.^2 .* ((1 - theta).^2) .* theta; % Función de utilidad

%% Paso 5: Configurar la iteración de la función de valor
dt = 0.1;   % Paso de tiempo fijo
tol = 1e-6; % Tolerancia de convergencia
max_iter = 1000; % Máximo de iteraciones

%% Definir theta inicial aleatorio
theta = rand(N, 1);

%% Iteración de la ecuación de Bellman
for iter = 1:max_iter
    % Matriz del método implícito con dt fijo
    I = eye(N);
    A = (1 + rho * dt) * I - dt * M; % Matriz a invertir

    % Evaluar u(theta, sigma) en cada iteración
    u_eval = u_theta(theta, x);
    
    % Resolver para theta en esta iteración
    theta_new = A \ (dt * u_eval + theta);
    
    % Criterio de convergencia
    diff = norm(theta_new - theta, Inf);
    if diff < tol
        fprintf('Convergencia alcanzada en %d iteraciones.\n', iter);
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

