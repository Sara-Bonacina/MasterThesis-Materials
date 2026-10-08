function p=ToxX1Dinit(p,lx,nx,par)
    %% setting standard parameters and screenlayout
    p=stanparam(p); % infuses p with standard parameter settings
    screenlayout(p); % open, clear and arrange the common figures
    %% special parameters related to this model
    % basics
    p.nc.neq=2; % number of equations in the model
    p.nc.nsteps=300;
    p.sw.bifcheck=2;
    p.sw.foldcheck=1;
    p.sw.sfem=-1; % type of numerical calculation, here OOPDE
    p.sw.jac=0; % use numerical jacobian (0) or analytical jacobian (1)
    p.sw.spjac=0; % use analytical Jacobian for spectral point cont (fold cont)
    p.plot.pcmp=1; % plotsol plots the 2nd component
    % description of the model
    p.fuha.sG=@sG; % the model itself
    %p.fuha.sGjac=@sGjac; % the Jacobian of the model  
    p.fuha.outfu=@ToxicityX1Dbra; % output quantities
    %% domain and mesh
    p.pdeo=stanpdeo1D(lx,2*lx/nx); % mesh [-lx ,lx], max mesh pt 2* lx/r
    bc=p.pdeo.grid.neumannBC('0');
    p.pdeo.grid.makeBoundaryMatrix(bc);
    p.np=p.pdeo.grid.nPoints ; % number of meshpoints
    p.nu=p.np*p.nc.neq; % number of unknowns (=2*( mesh points ), as 3 components )
    %%
    % construction the trivial solution
    % par=[g,gamma,Rc,d,s,c,k,dR,dT,psigma]';
    g=par(1); gamma=par(2); Rc=par(3); d=par(4); s=par(5); c=par(6); k=par(7);
    
    Tc=Rc*c*(d+s)/k; % sort of carrying capacity for T
    
    Rs = (g-d)*Rc/g;
    Ts = c*d*Rs/k;
        
    p.u=[Rs*ones(p.np,1);Ts*ones(p.np,1);par]; % initial solution guess with parameters
    %%
    p=setfemops(p); % compute FEM - operators
    %% bifurcation parameter , continuation basics and first guess for solution
    p.nc.ilam=5; % primary bifurcation parameter located at p.u(p.np+p.nc.ilam )
    p.nc.lammin=0; % lower bound for primary parameter during continuation
    p.nc.lammax=1; % upper bound for primary parameter during continuation
    p.usrlam =[0.5]; % user - vals for output
    p.sol.xi=1/p.nu; % weight in arclength - continuation
    p.sol.ds=-0.0001; % starting stepsize
    p.nc.dsmax=1e-4; % maximal stepsize
    p.nc.dsmin=1e-7; % minimal stepsize
    %% Plot
    p.plot.auxdict ={'g','\gamma','\hat{R}','d','s','c','k','d_R','d_T','\sigma/d_R'}; % parameter names 
end