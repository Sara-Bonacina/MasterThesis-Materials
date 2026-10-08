
clc
keep pphome 
%% bifurcation diagram - parameter sigma
% parameters kept fixed
g=10; % growth rate
d=1; % mortality rate                       ! g>d !
c=0.5; % conversion factor
k=1; % toxixity depletion
Rc=6; % sort of carrying capacity for R
dR=3.33; % diffusion parameter (biomass)
dT=0.05; % diffusion parameter (toxicity)

%% parameters varied
gamma= 0.1; % growth inhibition
s= 0; % extramortality
psigma_ref =0.9; % 

delta=0.1;
psigma = psigma_ref - delta; % starting point for continuation

par=[g,gamma,Rc,d,s,c,k,dR,dT,psigma]';

lx=4;  % domain (interval)
p=[];
p=ToxX1Dinit(p,lx,400,par); 
p.Om=2*lx;
p = setfn(p, 'cap3Case4'); 

p.nc.ilam=10; % continuation wrt p_sigma
%% continuation of the trivial branch
p.nc.dsmax=1e-3; 
p.nc.dsmin= 1e-6;   
p.sw.bifcheck=2;
p.sw.spcalc = 1;
p.sol.ds=0.01; % starting stepsize

p.nc.lammax=1; % sigma<dR 
p.nc.lammin=psigma; 
p=cont(p,1500); % continuation of the homogeneous branch

% File list
currentFolder=pwd;
cd('cap3Case4'); 

ptFileList = dir('bpt*.mat'); 

[nn,~]=size(ptFileList);
cd(currentFolder)
%%
sig_bif = zeros(nn,1);
for i = 1:nn
    S = load(fullfile('cap3Case4', ['bpt' num2str(i) '.mat']));
    q = S.p;
    sig_bif(i) = q.u(q.nu + q.nc.ilam(1));    % parameter index 10 -> psigma = sigma/d_R
end
[sig_sorted, idx] = sort(sig_bif);
fprintf('first bifurcation: sigma/d_R = %.4f  (sigma = %.4f)\n', sig_sorted(1), sig_sorted(1)*dR);
disp(table(idx, sig_sorted, sig_sorted*dR, 'VariableNames', {'bpt','sigma_over_dR','sigma'}));
%% continuation from bif points
for i=1:nn
    BPT_list=['bpt' num2str(i)];
    %Branch_i_u=['bpt' num2str(i) '_u'];
    Branch_i_d=['bpt' num2str(i) '_d'];
    p=swibra('cap3Case4',BPT_list,Branch_i_d,1e-2);  p.nc.dsmax=1e-1; p.sw.bifcheck=2; p.plot.pcmp=1; p=cont(p,100); p.nc.lammax=1;  
    %p=swibra('prova',BPT_list,Branch_i_u,1e-2);  p.nc.dsmax=1e-2; p.sw.bifcheck=1; p.plot.pcmp=1; p=cont(p,500); p.nc.lammax=1;  
    
end

%% Postprocessing, plot BifDiagram
nfig=5;
figure(nfig);
clf(nfig);
cmp=2;
box on
hold on

col_main = [0 0 0];              % black - homog. branch
col_sec  = [0 0.4471 0.6980];    % blue - secondary branches
% homog. branch in black (thick stable, thin unstable) 
plotbra('cap3_case4',nfig,cmp,'cl','col_main','tyst','-','tyun','--');

% secondary branches in blue(thick stable/ thin unstable)
for i=[1, 2:nn]
    Branch_i_d=['bpt' num2str(i) '_d'];
    plotbra(Branch_i_d,nfig,cmp,'cl',col_sec, 'tyst','-','tyun','-');
end

xlabel('\sigma/d_R')
ylabel('||R||_{L^1}')
%axis([0 0.85 302 347])

%% Select branch for video
branch_idx = 3;      % <-- indice della biforcazione di interesse
branch_dir = 'd';    % <-- 'u' o 'd'
col_sel  = [0.8353 0.3686 0];    % vermiglio - ramo selezionato


branchOfInterest = ['bpt' num2str(branch_idx) '_' branch_dir];
plotbra(branchOfInterest,nfig,cmp,'cl',col_sel, 'tyst','-','tyun','--');

savefig(figure(nfig), 'BD3_psigma_case4.fig');
xlabel('\sigma/d_R')
ylabel('||R||_{L^1}')

% copy the folder of the branch for SolMovie
targetFolder = 'Branch3_case4';
if exist(targetFolder,'dir')
    rmdir(targetFolder,'s');
end
copyfile(branchOfInterest, targetFolder);

fprintf('Branch "%s" salvato in "%s" (%d file .mat)\n', ...
    branchOfInterest, targetFolder, numel(dir(fullfile(targetFolder,'pt*.mat'))));


