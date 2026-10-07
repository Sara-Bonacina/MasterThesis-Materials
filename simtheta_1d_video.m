% 1d simulation with different approximation of theta (to select at the
% end, in the pde function) with temporal evolution and video saved (change
% the name)
% comparison of 4 scenarios with effects of s, sigma, gamma
% power function: no patterns in the last case, try with a different I.C/ different sigma (e.g. sigma=0.95dR should works)!!!

clear all; close all; clc;
tic;

% space discretization
L = 8; dx = 0.02; x = 0:dx:L; Nx = length(x);
% time discretization
t_start = 0; t_end = 100;
t_video = linspace(t_start, t_end, 200); % 200 frame per il video
tspan = [t_start t_end];

% Parameters of the model: 4 cases to see influence of s, gamma, sigma
g = 10; c = 0.5; Rc = 6; d = 1; k = 1;
diff_R = 3.33; diff_T = 0.05;
opt = [0.1, 0.5, 0; 0.1, 0.5, 0.9*diff_R; 0, 0.5, 0.9*diff_R; 0.1, 0, 0.9*diff_R];

% epsilon smooth approx
epsilon = 0.02;

% power function exponents
a = 5; b = 4;

% for solutions
sol_R = cell(4,1);
sol_T = cell(4,1);

options = odeset('RelTol', 1e-4, 'AbsTol', 1e-6);
for index = 1:4
    p.gamma = opt(index,1); p.s = opt(index,2); p.sigma = opt(index,3);
    p.Tc = c*(d + p.s)*Rc/k;
    p.g = g; p.c = c; p.Rc = Rc; p.d = d; p.k = k; p.dR = diff_R; p.dT = diff_T;
    p.epsilon = epsilon;
    p.a = a; p.b = b;
    p.dx = dx; p.Nx = Nx;
    
    % % I.C: try also with a different i.c (e.g. perturbation of equilibrium)
    % R0 = 10*exp((-1/L)*(5*x - 5*L/2).^2)'; 
    % T0 = zeros(Nx, 1);

    % % I.C : equilibrium perturbation
    R0 = ((g-d)*Rc/g)*(1+0.8*(2*rand(Nx,1)-1));
    T0 = (c*d*R0/k);
    y0 = [R0; T0];
 
    % 
    [t_out, sol_out] = ode15s(@(t, y) system_pde(t, y, p), tspan, [R0; T0], options);
    
    % Interpolation on temporal grid of the video for uniformity
    % sol_out has dimension [Time x (2*Nx)]
    sol_interp = interp1(t_out, sol_out, t_video);
    
    sol_R{index} = sol_interp(:, 1:Nx);
    sol_T{index} = sol_interp(:, Nx+1:end);
    fprintf('Case %d completed\n', index);
end

% 2 Video (GIF): DON'T FORGET TO CHANGE THE NAME!
nomeFile = 'evolution_ODE15s.gif';
hFig = figure('Position', [100 100 1000 800], 'Color', 'w');
y_max = 10;

for k_frame = 1:length(t_video)
    for index = 1:4
        subplot(2,2,index);
        % data of the current frame for the specific case
        R_curr = sol_R{index}(k_frame, :);
        T_curr = sol_T{index}(k_frame, :);
        
        plot(x, R_curr, 'b-', 'LineWidth', 1.5); hold on;
        plot(x, T_curr, 'r--', 'LineWidth', 1.5); hold off;
        
        ylim([0 y_max]); grid on;
        title(sprintf('Case %d: \\sigma=%.2f (t=%.1f)', index, opt(index,3), t_video(k_frame)));
        if index == 1, legend('R','T'); end
    end
    
    drawnow;
    
    % Save the frame into the GIF
    frame = getframe(hFig);
    im = frame2im(frame);
    [imind, cm] = rgb2ind(im, 256);
    if k_frame == 1
        imwrite(imind, cm, nomeFile, 'gif', 'Loopcount', inf, 'DelayTime', 0.05);
    else
        imwrite(imind, cm, nomeFile, 'gif', 'WriteMode', 'append', 'DelayTime', 0.05);
    end
end

toc;

%% PDE function
function dydt = system_pde(t, y, p)
    Nx = p.Nx; dx = p.dx;
    Tc = p.Tc;
    epsilon = p.epsilon;
    a=p.a; b=p.b;
    R = y(1:Nx);
    T = y(Nx+1:end);
    
    % DON'T FORGET TO SELECT THETA!:
    %%% theta original
    theta = T / p.Tc;
    theta(theta > 1) = 1;
    %%% smooth approx:
    % theta = 1-epsilon*log(1+exp((Tc-T)./(Tc*epsilon)));
    %%% power function 
    % Ts = T/Tc;
    % theta = (Ts.^a) ./ (Ts.^a + (1+Ts).^(-b));

    dRdt = zeros(Nx, 1);
    dTdt = zeros(Nx, 1);

    % cross diffusion
    ArgR = (p.dR - p.sigma * theta) .* R;
    
    % Neumann BC idx=1
    idx = 1;
    lapR = (ArgR(idx+1) - 2*ArgR(idx) + ArgR(idx+1)) / p.dx^2;
    lapT = p.dT * (T(idx+1) - 2*T(idx) + T(idx+1)) / p.dx^2;
    
    % Reaction
    reacR = (p.g - p.gamma*theta(idx)) .* R(idx) .* (1 - R(idx)/p.Rc) ...
            - (p.d + p.s*theta(idx)) .* R(idx);
    reacT = p.c * (p.d + p.s*theta(idx)) .* R(idx) - p.k * T(idx);
    
    dRdt(idx) = lapR + reacR;
    dTdt(idx) = lapT + reacT;

    % Internal index
    i = 2:Nx-1;
    
    % Cross-diffusion term: diff( (dR - sigma*theta)*R )
    argdiff = (p.dR - p.sigma * theta) .* R;
    diffR = (argdiff(i+1) - 2*argdiff(i) + argdiff(i-1)) / dx^2;
    diffT = p.dT * (T(i+1) - 2*T(i) + T(i-1)) / dx^2;
    
    % Reaction terms
    reacR = p.g*R(i).*(1 - R(i)/p.Rc) - ...
            (p.gamma*theta(i)).*R(i).*(1 - R(i)/p.Rc) - ...
            (p.d + p.s*theta(i)).*R(i);
            
    reacT = p.c*(p.d + p.s*theta(i)).*R(i) - p.k*T(i);
    
    dRdt(i) = reacR + diffR;
    dTdt(i) = reacT + diffT;
    
    % Neumann BC idx=Nx
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