% figure 2: 1d simulation 4 scenarios with ode15s solver. 
% NO temporal evolution, 30 seconds to plot:
% I.C., 4 scenarios B and T, 1 plot for B, T in the 4 cases at final time (t=1000). 

clear all; close all; clc;
%
tic;

% space and time
L = 8;
dx = 0.02;
x = 0:dx:L;
Nx = length(x);
tspan = [0 1000];

% parameters
g = 10; c = 0.5; Rc = 6; d = 1; k = 1;
diff_R = 3.33; diff_T = 0.05;

% (gamma, s, sigma)
opt = [0.1, 0.5,   0; 
       0.1,   0.5, 0.9*diff_R; 
       0, 0.5, 0.9*diff_R; 
       0.1, 0, 0.9*diff_R];

% final results
R_all = zeros(size(opt,1), Nx);
T_all = zeros(size(opt,1), Nx);

figure('Position', [100 100 1200 800]);

y_max = 10;

for index = 1:size(opt,1)
    params.gamma = opt(index,1);
    params.s     = opt(index,2);
    params.sigma = opt(index,3);
    params.Tc    = c*(d + params.s)*Rc/k;
    params.g = g; params.c = c; params.Rc = Rc;
    params.d = d; params.k = k; params.dR = diff_R; params.dT = diff_T;
    params.dx = dx; params.Nx = Nx;

    % I.C.
    R0 = 10*exp((-1/L)*(5*x - 5*L/2).^2)'; 
    T0 = zeros(Nx, 1);
    y0 = [R0; T0]; 

    % ode15s 
    options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
    [t, sol] = ode15s(@(t, y) system_pde(t, y, params), tspan, y0, options);

    % 
    R_final = sol(end, 1:Nx);
    T_final = sol(end, Nx+1:end);

    % Plot
    subplot(2,2,index);
    plot(x, R_final, 'b-', 'LineWidth', 1.5); 
    hold on;
    plot(x, T_final, 'r--', 'LineWidth', 1.5); 
    hold off;

    % fix limits for y-axes
    ylim([0 y_max]); 
    ylabel('R , T');
    title(sprintf('Case %d: \\gamma=%.2f, s=%.2f, \\sigma=%.2f ', index, params.gamma, params.s, params.sigma));
    grid on;

    R_all(index, :) = sol(end, 1:Nx);
    T_all(index, :) = sol(end, Nx+1:end);
    fprintf('end case %d\n', index);    
end

% plot I.C.:
figure('Name', 'Initial condition', 'Position', [100, 100, 1200, 500]);
plot(x, R0, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Biomass (R)');
hold on;
plot(x, T0, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Toxicity (T)');
hold off;
title('Initial condition');
xlabel('x');
ylabel('R, T');
legend('Location', 'northeast');
grid on;

% Plot to compare 
figure('Name', 'Comparison 4 Cases', 'Position', [100, 100, 1200, 500]);
nomi_casi = {'(i): \gamma, s', '(ii): \gamma, s, \sigma', '(iii): s, \sigma', '(iv): \gamma, \sigma'};
colori = {'b', 'r', 'g', 'k'}; % different colors for 4 cases
stili = {'-', '--', ':', '-.'}; % different style for lines

% Subplot 1: Biomass R
subplot(1, 2, 1);
hold on;
for i = 1:4
    plot(x, R_all(i,:), 'Color', colori{i}, 'LineStyle', stili{i}, 'LineWidth', 1.5);
end
title('Biomass R');
xlabel('x'); ylabel('R');
legend(nomi_casi, 'Location', 'best');
grid on;

% Subplot 2: Toxicity T
subplot(1, 2, 2);
hold on;
for i = 1:4
    plot(x, T_all(i,:), 'Color', colori{i}, 'LineStyle', stili{i}, 'LineWidth', 1.5);
end
title('Toxicity T');
xlabel('x'); ylabel('T');
legend(nomi_casi, 'Location', 'best');
grid on;

hold off;
toc

%--------------------------------------------------------------------------
% PDE discretized
function dydt = system_pde(~, y, p)
    Nx = p.Nx;
    
    R = y(1:Nx);
    T = y(Nx+1:end);
   
    dRdt = zeros(Nx, 1);
    dTdt = zeros(Nx, 1);
    
    % theta, argR for crosss-diff
    theta = min(T / p.Tc, 1);
    ArgR = (p.dR - p.sigma * theta) .* R;
    
    % Neumann B.C. idx=1
    idx = 1;
    lapR = (ArgR(idx+1) - 2*ArgR(idx) + ArgR(idx+1)) / p.dx^2;
    lapT = p.dT * (T(idx+1) - 2*T(idx) + T(idx+1)) / p.dx^2;
    
    % Reaction
    reacR = (p.g - p.gamma*theta(idx)) .* R(idx) .* (1 - R(idx)/p.Rc) ...
            - (p.d + p.s*theta(idx)) .* R(idx);
    reacT = p.c * (p.d + p.s*theta(idx)) .* R(idx) - p.k * T(idx);
    
    dRdt(idx) = lapR + reacR;
    dTdt(idx) = lapT + reacT;
    
    % internal nodes
    idx = 2:Nx-1;
    
    % diffusion
    lapR = (ArgR(idx+1) - 2*ArgR(idx) + ArgR(idx-1)) / p.dx^2;
    lapT = p.dT * (T(idx+1) - 2*T(idx) + T(idx-1)) / p.dx^2;
    
    % Reaction
    reacR = (p.g - p.gamma*theta(idx)) .* R(idx) .* (1 - R(idx)/p.Rc) ...
            - (p.d + p.s*theta(idx)) .* R(idx);
    reacT = p.c * (p.d + p.s*theta(idx)) .* R(idx) - p.k * T(idx);
    
    dRdt(idx) = lapR + reacR;
    dTdt(idx) = lapT + reacT;

    % Neumann B.C. idx=Nx
    idx = Nx;
    lapR = (ArgR(Nx-1) - 2*ArgR(idx) + ArgR(idx-1)) / p.dx^2;
    lapT = p.dT * (T(Nx-1) - 2*T(idx) + T(idx-1)) / p.dx^2;
    
    % Reaction
    reacR = (p.g - p.gamma*theta(idx)) .* R(idx) .* (1 - R(idx)/p.Rc) ...
            - (p.d + p.s*theta(idx)) .* R(idx);
    reacT = p.c * (p.d + p.s*theta(idx)) .* R(idx) - p.k * T(idx);
    
    dRdt(idx) = lapR + reacR;
    dTdt(idx) = lapT + reacT;
      
    dydt = [dRdt; dTdt];
end