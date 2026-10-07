% simulation 2D in the 4 scenarios
% L=8, lattice 80x80, dx=dy=0.1, N.B.C.
% I.C.: 4 squares of 4 pixels per side randomly placed in the lattice
% with uniform value of R randomly chosen in range 1-2 kg/m^2, T=0
%
% Saves:
%   1. Initial condition
%   2. Final state
%   3. 100 frames for Beamer animation

clear all; close all; clc;
tic

% parameters
g = 10; 
c = 0.5; 
Rc = 6; 
d = 1; 
k = 1;

diff_R = 3.33; 
diff_T = 0.05;

% Spatial discretization
L = 8;
dx = 0.1; 
dy = 0.1;
x = 0:dx:L; 
y = 0:dy:L;
Nx = length(x); 
Ny = length(y);

% Time discretization
dt = 0.0005;
t_end = 300;

% 
sigma_on = 0.9*diff_R;

cases = struct( ...
    'label', {'(i) No propagation-reduction', ...
              '(ii) All effects active', ...
              '(iii) No growth-inhibition', ...
              '(iv) No extra-mortality'}, ...
    'gamma', {0.1, 0.1, 0,   0.1}, ...
    's',     {0.5, 0.5, 0.5, 0}, ...
    'sigma', {0,   sigma_on, sigma_on, sigma_on} ...
);

nCases = numel(cases);

% I C
rng(1);

R0 = zeros(Ny,Nx);

side = 4;

for q = 1:4
    ix = randi([2,Nx-side]);
    iy = randi([2,Ny-side]);

    val_random = 1 + rand();

    R0(iy:iy+side-1,ix:ix+side-1) = val_random;
end

T0 = zeros(Ny,Nx);

R = cell(1,nCases);
T = cell(1,nCases);

R_final_all = cell(1,nCases);
T_final_all = cell(1,nCases);

for cidx = 1:nCases
    R{cidx} = R0;
    T{cidx} = T0;
end

% Video settings
video_folder = 'frames_beamer_2D_4cases';

if ~exist(video_folder,'dir')
    mkdir(video_folder);
end

% 100 frames, with more resolution during the initial evolution
t1 = linspace(0,50,42);
t2 = linspace(50,150,35);
t3 = linspace(150,300,27);

t_video = [t1,t2(2:end),t3(2:end)];

n_frames = length(t_video);
next_frame = 1;

% Figure for animation
FigVideo = figure( ...
    'Position',[50 50 1400 700], ...
    'Color','w');

% loop
Time = 0;
for t_step = 0:dt:t_end
    for cidx = 1:nCases

        gamma = cases(cidx).gamma;
        s     = cases(cidx).s;
        sigma = cases(cidx).sigma;

        Tc = c*(d+s)*Rc/k;

        R_old = R{cidx};
        T_old = T{cidx};

        theta = T_old/Tc;
        theta(theta > 1) = 1;

        %Cross-diffusion
        argdiff = (diff_R - sigma*theta).*R_old;
        % Neumann BC
        argdiffPad = zeros(Ny+2,Nx+2);

        argdiffPad(2:end-1,2:end-1) = argdiff;

        argdiffPad(1,2:end-1)   = argdiff(2,:);
        argdiffPad(end,2:end-1) = argdiff(end-1,:);
        argdiffPad(2:end-1,1)   = argdiff(:,2);
        argdiffPad(2:end-1,end) = argdiff(:,end-1);

        TPad = zeros(Ny+2,Nx+2);

        TPad(2:end-1,2:end-1) = T_old;

        TPad(1,2:end-1)   = T_old(2,:);
        TPad(end,2:end-1) = T_old(end-1,:);
        TPad(2:end-1,1)   = T_old(:,2);
        TPad(2:end-1,end) = T_old(:,end-1);

        % Laplacian of R
        LapR = ...
            (argdiffPad(2:end-1,3:end) ...
            - 2*argdiffPad(2:end-1,2:end-1) ...
            + argdiffPad(2:end-1,1:end-2))/dx^2 ...
            + ...
            (argdiffPad(3:end,2:end-1) ...
            - 2*argdiffPad(2:end-1,2:end-1) ...
            + argdiffPad(1:end-2,2:end-1))/dy^2;

        % Laplacian of T
        LapT = diff_T*( ...
            (TPad(2:end-1,3:end) ...
            - 2*TPad(2:end-1,2:end-1) ...
            + TPad(2:end-1,1:end-2))/dx^2 ...
            + ...
            (TPad(3:end,2:end-1) ...
            - 2*TPad(2:end-1,2:end-1) ...
            + TPad(1:end-2,2:end-1))/dy^2 );

        % Reaction terms
        reacR = ...
            g*R_old.*(1-R_old/Rc) ...
            - (gamma*theta).*R_old.*(1-R_old/Rc) ...
            - (d+s*theta).*R_old;

        reacT = ...
            c*(d+s*theta).*R_old-k*T_old;

       % forward euler
        R{cidx} = R_old + (reacR+LapR)*dt;
        T{cidx} = T_old + (reacT+LapT)*dt;
    end

    Time = Time + dt;
  
    % saving frame
    if next_frame <= n_frames && Time >= t_video(next_frame)
        figure(FigVideo);
        for cidx = 1:nCases
            % biomass
            ax = subplot(2,4,cidx);
            imagesc(x,y,R{cidx});
            axis square;
            set(gca,'YDir','normal');
            caxis([0 10]);
            title(cases(cidx).label, ...
                'FontSize',10, ...
                'Interpreter','latex');
            set(gca,'xtick',[],'ytick',[]);
            colormap(ax,flipud(summer));

            if cidx == 1
                ylabel('Biomass $R$', ...
                    'FontSize',12);
            end

            colorbar;
            % toxicity
            ax = subplot(2,4,4+cidx);
            imagesc(x,y,T{cidx});
            axis square;
            set(gca,'YDir','normal');
            caxis([0 5]);
            set(gca,'xtick',[],'ytick',[]);
            colormap(ax,flipud(gray));

            if cidx == 1
                ylabel('Toxicity $T$', ...
                    'FontSize',12);
            end

            colorbar;

        end

        sgtitle(sprintf('$t = %.1f$',Time), ...
            'FontSize',14);

        drawnow;

        filename = fullfile( ...
            video_folder, ...
            sprintf('frame_%03d.png',next_frame));

        exportgraphics(FigVideo, ...
            filename, ...
            'Resolution',150);

        fprintf('Frame %d/%d saved: t = %.2f\n', ...
            next_frame,n_frames,Time);

        next_frame = next_frame + 1;

    end

end

% Save final states

for cidx = 1:nCases
    R_final_all{cidx} = R{cidx};
    T_final_all{cidx} = T{cidx};
end

close(FigVideo);

toc


%  FIGURE 1: Initial condition
Fig1 = figure('Position',[50 50 700 350]);

ax = subplot(1,2,1);

imagesc(x,y,R0);
axis square;
set(gca,'YDir','normal');

caxis([0 10]);

title('Biomass $R$', ...
    'FontSize',12);

set(gca,'xtick',[],'ytick',[]);

colormap(ax,flipud(summer));
colorbar;
ax = subplot(1,2,2);

imagesc(x,y,T0);
axis square;
set(gca,'YDir','normal');

caxis([0 5]);

title('Toxicity $T$', ...
    'FontSize',12);

set(gca,'xtick',[],'ytick',[]);

colormap(ax,flipud(gray));
colorbar;


sgtitle('Initial condition ($t=0$)', ...
    'FontSize',14);


% FIGURE 2: Final state

Fig2 = figure('Position',[50 50 1400 700]);

for cidx = 1:nCases
    % Biomass
    
    ax = subplot(2,nCases,cidx);

    imagesc(x,y,R_final_all{cidx});

    axis square;
    set(gca,'YDir','normal');

    caxis([0 10]);

    title(cases(cidx).label, ...
        'FontSize',10);

    set(gca,'xtick',[],'ytick',[]);

    colormap(ax,flipud(summer));

    if cidx == 1
        ylabel('Biomass $R$', ...
            'FontSize',12);
    end

    if cidx == nCases
        pos = get(ax,'Position');
        colorbar;
        set(ax,'Position',pos);
    end

    % Toxicity
   
    ax = subplot(2,nCases,nCases+cidx);

    imagesc(x,y,T_final_all{cidx});

    axis square;
    set(gca,'YDir','normal');

    caxis([0 5]);

    set(gca,'xtick',[],'ytick',[]);

    colormap(ax,flipud(gray));

    if cidx == 1
        ylabel('Toxicity $T$', ...
            'FontSize',12);
    end

    if cidx == nCases
        pos = get(ax,'Position');
        colorbar;
        set(ax,'Position',pos);
    end

end

sgtitle(sprintf('Final state ($t=%d$)',t_end), ...
    'FontSize',14);
