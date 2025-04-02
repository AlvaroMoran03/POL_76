function solve_ode
    % Parámetros
    sigma = 0.1;   % Define un valor adecuado para sigma
    rho_e = 0.05;  % Parámetro de la ecuación
    rho_h = 0.03;  % Parámetro de la ecuación
    q_h = 1;       % Condición de frontera
    eta_span = [0 1]; % Intervalo de integración

    % Resolver la ODE con la condición inicial q(0) = q_h
    [eta, q] = ode45(@(eta, q) ode_system(eta, q, sigma, rho_e, rho_h), eta_span, q_h);
    
    % Graficar la solución
    figure;
    plot(eta, q, 'b', 'LineWidth', 2);
    xlabel('\eta');
    ylabel('q(\eta)');
    title('Solución de la ODE');
    grid on;
end

function dq_deta = ode_system(eta, q, sigma, rho_e, rho_h)
    % Definición de kappa^e
    a_e = 1; % Define un valor adecuado
    a_h = 1; % Define un valor adecuado
    l_q = q; % Suponemos l(q) = q para simplificar
    
    kappa_e = (l_q + q * (eta * rho_e + (1 - eta) * rho_h) - a_h) / (a_e - a_h);
    
    % Definir la ecuación diferencial
    if kappa_e < 1
        dq_deta = sqrt((sigma / (1 - (kappa_e / eta - 1)))^2); % Se toma la raíz positiva
    else
        dq_deta = 0; % Condición cuando kappa^e = 1
    end
end
