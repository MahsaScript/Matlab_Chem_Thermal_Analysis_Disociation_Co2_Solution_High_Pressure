% This program may be copied and used in unaltered form. Changes may be
% made to the contents at any point for personal use. Redistribution of
% this program in both altered and unaltered form is subject to the
% conditions that it remains free of charge and that this notice is
% included upfront. No warranties are given, expressed or implied, that
% this program is free of errors and bugs. Neither do we guarantee that
% this program meets the requirements of your particular application. The
% author disclaims all liability for direct or consequential damages as a
% result from use of this program.

function tdr_model(T,P,j,n,electrolyte,c,L,n_x,t,dt,simplified)
% This model uses the pdepe.m solver to find numerical solutions to the
% governing transport equations of the form
%
%   d/dt(c) = D*d^2/dx^2(c) + R_V 
%
% which describes the transient reaction-diffusion equation for species x
% with volumetric reaction terms R_V (which couples the various equations).
% This is a 1D system describing the diffusion layer between x = 0 and x = 
% L, with the electrode positioned at x = 0 and the bulk solution at x = L.
%
% Note that the concentrations returned by the solver are in units of
% [mol/m3] (i.e. [mM]) to keep in line with the dimensionality of the
% governing equations. For equilibrium and rate coefficients, their values
% are outputted in the console window in dimensions involving [mol/L] for
% ease of comparison with literature values, after which they are converted
% to [mol/m3]-equivalent values and used as such in the solution procedure.

    clc; close all;
    format compact; format short;
    fprintf('***1D Transient Diffusion-Reaction Model for CO2 Electrolysis***\n\n')

    % cast input parameters in other units
    j             = j*10;                             % current density [A/m2]
    T             = 273.15 + T;                       % temperature [K]
    x_0           = 0;                                % model domain left boundary [m]
    x_DBL         = L/1e6;                            % diffusion layer thickness [m]
    % values of other natural constants
    F             = 96485;                            % Faraday's constant [C/mol]
    % plot colours
    b = [0 0.4470 0.7410];
    o = [0.8500 0.3250 0.0980];
    y = [0.9290 0.6940 0.1250];
    g = [0.4660 0.6740 0.1880];

    tic
    % (1) Initialize parameters relevant to the model.
    fprintf('\n---PART I: Thermodynamic and Kinetic Parameters---\n\n');
    % viscosity and density of KHCO3 electrolyte in [Pa s] as function of molality
    [u,p]         = viscodensi(T,c/1e3,electrolyte);
    fprintf('Thermodynamic state of the system:\n');
    fprintf('\t T = %.2f K\n',T);
    fprintf('\t P = %.2f bar\n',P);
    fprintf('\t c = %.2f mol/m3 %s\n',[c,electrolyte]);
    fprintf('\t u = %.2f mPa s (based on salt molarity and temperature)\n',u*1e3);
    fprintf('\t p = %.1f kg/m3 (based on salt molarity)\n',p);
    % saturated CO2 concentration [mol/m3] based on electrolyte and T,P conditions 
    c_sat = sechenov(c/1e3,T,P,'CO2',electrolyte)*1e3;
    fprintf('Solubility of CO2 in equilibrium with %.2f bar pure CO2 atmosphere:\n',P);
    fprintf('\t %.1f mM (Henry Law)\t --->\t %.1f mM (Sechenov)\n',henry(T,P,'CO2')*1e3,c_sat);

    % diffusivity of species in aqueous phase
    D_DBL = diffusivity(T,u,0,1);
    fprintf('Diffusion coefficients of species in the diffusion layer:\n')
    fprintf('\t D(CO2)     = %2.3e m2/s\n',D_DBL(1).CO2);
    fprintf('\t D(HCO3(-)) = %2.3e m2/s\n',D_DBL(1).iHCO3);
    fprintf('\t D(CO3(2-)) = %2.3e m2/s\n',D_DBL(1).iCO3);        
    fprintf('\t D(OH(-))   = %2.3e m2/s\n',D_DBL(1).iOH);
    fprintf('\t D(H(+))    = %2.3e m2/s\n',D_DBL(1).iH);
    fprintf('Corrected for temperature and viscosity.\n\n');
    
    fprintf('Self-ionization constant of water:\n');
    pKw = selfionization(T,0,electrolyte);
    Kw  = 10^(-pKw);       
    fprintf('\t Kw         = %2.3e M2 (pKw = -log(Kw) = %1.2f)\n',[Kw, pKw]);
    fprintf('Note that Kw is only corrected for temperature, not ionic strength.\n\n');

    fprintf('Equilibrium coefficients:\n');
    %   CO2 + H2O <---> HCO3(-) + H(+)  with Ka1 in [mol/L]
    %   HCO3(-)   <---> CO3(2-) + H(+)  with Ka2 in [mol/L]
    %   HCO3(-)   <---> CO2 + OH(-)     with Kb1 in [mol/L]
    %   CO3(2-)   <---> HCO3(-) + OH(-) with Kb2 in [mol/L]
    [Ka1, Ka2, ~, ~] = carbonateeq(T,c/1e3,electrolyte);
    Kb1              = Kw/Ka1;
    Kb2              = Kw/Ka2;
    fprintf('\t Ka1        = [HCO3(-)][H(+)] /[CO2]     = %2.3e M (pKa1 = %1.2f)\n',[Ka1, -log10(Ka1)]);
    fprintf('\t Ka2        = [CO3(2-)][H(+)] /[HCO3(-)] = %2.3e M (pKa2 = %1.2f)\n',[Ka2, -log10(Ka2)]);
    fprintf('\t Kb1        =     [CO2][OH(-)]/[HCO3(-)] = %2.3e M (pKb1 = %1.2f)\n',[Kb1, -log10(Kb1)]);
    fprintf('\t Kb2        = [HCO3(-)][OH(-)]/[CO3(2-)] = %2.3e M (pKb2 = %1.2f)\n',[Kb2, -log10(Kb2)]);
    fprintf('Corrected for temperature and ionic strength.\n\n');

    % Sources report k-values involving [L] units. These are reported here
    % for ease of comparison, before being converted to units involving [m]
    % to keep in line with the dimensions of the governing equations.
    [k3,k5] = carbonatekin(T,c/1e3,electrolyte);
    k       = struct('ka1f', {k3},...
                     'ka1b', {k3/Ka1},...
                     'ka2f', {Ka2*5.0e10},...
                     'ka2b', {5.0e10},...
                     'kb1f', {k5/Kw*Kb1},...
                     'kb1b', {k5/Kw},...
                     'kb2f', {6e9*Kb2},...
                     'kb2b', {6e9},...
                     'kwf', {1.4e-3},...
                     'kwb', {1.4e-3/Kw});
    fprintf('Rate coefficients:\n');
    fprintf('\t ka1>        = %2.3e /s\t\t ka1<        = %2.3e /M/s\n',[k.ka1f, k.ka1b]);
    fprintf('\t ka2>        = %2.3e /s\t\t ka2<        = %2.3e /M/s\n',[k.ka2f, k.ka2b]);
    fprintf('\t kb1>        = %2.3e /s\t\t kb1<        = %2.3e /M/s\n',[k.kb1f, k.kb1b]);
    fprintf('\t kb2>        = %2.3e /s\t\t kb2<        = %2.3e /M/s\n',[k.kb2f, k.kb2b]);
    fprintf('\t kw>        = %2.3e M/s\t\t kw<        = %2.3e /M/s\n',[k.kwf, k.kwb]);
    fprintf('Corrected for temperature and ionic strength.\n\n');
    % rewrite all kinetic coefficients in units involving [m3]
    k(1).ka1b = k(1).ka1b/1e3; k(1).ka2b = k(1).ka2b/1e3; k(1).kb1b = k(1).kb1b/1e3;
    k(1).kb2b = k(1).kb2b/1e3; k(1).kwf  = k(1).kwf*1e3;  k(1).kwb  = k(1).kwb/1e3;
    t1 = toc; fprintf('Elapsed time for PART I is %1.3f seconds.\n\n',t1);

    tic
    % (2) Calculation of concentrations in the bulk solution.
    fprintf('\n---PART II: Bulk Electrolyte Concentrations---\n');
    figure(1);
    [~,~,~,~] = bjerrum(Kw,Kb1,Kb2);
    title('Bjerrum plot of the carbonate system');
    fprintf('\nFraction of CO2/HCO3(-)/CO3(2-) of Dissolved Inorganic Content (DIC) as function of pH.\n');
    % using a root finding algorithm to find the initial conditions
    p_coeff   = [2*c_sat/1e3 (1 + c_sat/1e3/Kb1)*Kb1*Kb2 -c/1e3*Kb1*Kb2 -Kw*Kb1*Kb2];
    r         = roots(p_coeff);
    c_ic      = struct('CO2',   {c_sat},...
                       'iHCO3', {c_sat*r(2)/Kb1},...
                       'iCO3',  {c_sat*r(2)/Kb1*r(2)/Kb2},...
                       'iOH',   {r(2)*1e3},...
                       'iH',    {Kw/r(2)*1e3});
    pH_i      = -log10(c_ic.iH/1e3);
    fprintf('\nThe initial conditions of the CO2-saturated electrolyte are:\n');
    fprintf('\tCO2     %.5f mM\n',c_sat);
    fprintf('\tOH(-)   %.5f mM\n',c_ic(1).iOH);
    fprintf('\tHCO3(-) %.5f mM\n',c_ic(1).iHCO3);
    fprintf('\tCO3(2-) %.5f mM\n',c_ic(1).iCO3);
    fprintf('\tH(+)    %.5f mM\n',c_ic(1).iH);
    fprintf('\tpH      %.2f \n',pH_i);
    t2 = toc; fprintf('\nElapsed time for PART II is %1.3f seconds.\n\n',t2);

    tic
    % (3) Setting up and solving the system of matrix-vector equations.

    % set up logarithmic grid spacing
    dx           = -1/j*c_ic(1).iOH*F*D_DBL(1).iOH;
    x_left       = horzcat(x_0,logspace(log10(dx),log10(x_DBL/2),(n_x - 1)/2));
    x_right      = x_DBL - x_left; x_right = x_right(1:end - 1); x_right = flip(x_right);
    X            = horzcat(x_left,x_right);
    N            = [0:dt:t];

    fprintf('\n---PART III: Transient and Stabilized Solution Profiles---\n');
    fprintf('\nModel parameters:\n');
    fprintf('\tdiffusion boundary layer thickness     %.2f um\n',x_DBL*1e6);
    fprintf('\tlogarithmic grid spacing               %.0f points\n',n_x);
    fprintf('\tminimum step size at edge              %.5f um\n',dx*1e6)
    fprintf('\ttotal modelling time                   %.5f s\n',t);
    fprintf('\ttime step size                         %.5f s\n',dt);
    if simplified == 1
        fprintf('\tequations simplified\n');
    end
    
    % estimating the mass-transfer limiting current density.
    j_lim = -F*D_DBL(1).CO2*c_sat./x_DBL/(1/2*n(1) + 1/2*n(2) + 1/8*n(3) + 2/12*n(4) + 2/12*n(5) + 2/8*n(6));
    fprintf('\nFor j = %.1f mA/cm2 and n_CO = %.2f, the mass-transfer limiting current density is estimated as %1.1f mA/cm2.\n',[j/10,n(1),j_lim/10]);

    % solve the system of equations
    param_tdr    = @(x,t,u,dudx)tdr(x,t,u,dudx,D_DBL,k,simplified);
    param_tdr_ic = @(x)tdr_ic(x,c_ic);
    param_tdr_bc = @(xl,ul,xr,ur,t)tdr_bc(xl,ul,xr,ur,t,c_ic,j,F,n);
    sol     = pdepe(0,param_tdr,param_tdr_ic,param_tdr_bc,X,N);
    % dim 1: time T; dim 2: position X; dim 3: species type
    c_CO2   = sol(:,:,1);
    c_iHCO3 = sol(:,:,2);
    c_iCO3  = sol(:,:,3);
    c_iOH   = sol(:,:,4); 
    c_iH    = sol(:,:,5);
    pKw     = -log10(c_iH.*c_iOH/1e6);
    pH      = pKw + log10(c_iOH/1e3);

    % checking if the CO2 concentration is negative.
    if c_CO2(end,1) < 0
        fprintf('\nWarning: The partial current density exceeds the maximum supportable for CO2 reduction.\n')
    end

    t3 = toc; fprintf('\nElapsed time for PART III is %1.3f seconds.\n\n',t3);

    % (4) Preparing plots of the model solutions.
    tic

    figure(2);
    tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
    nexttile; plotting3D(X,N,c_CO2,'CO_{2}','{\it c} (mM)');
    nexttile; plotting3D(X,N,c_iHCO3,'HCO_{3}^{-}','{\it c} (mM)');
    nexttile; plotting3D(X,N,c_iCO3,'CO_{3}^{2-}','{\it c} (mM)');
    nexttile; plotting3D(X,N,c_iOH,'OH^{-}','{\it c} (mM)');
    nexttile; plotting3D(X,N,c_iH,'H^{+}','{\it c} (mM)');
    nexttile; plotting3D(X,N,pH,'pH','pH (-)');
    sgtitle('transient concentration profiles');
    fprintf('The transients upon imposing %1.1f mA/cm2 at the cathode.\n', j/10);

    % Plotting the stabilized solution profiles.
    figure(3);
    tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
    nexttile; plotting2D(X,c_CO2(end,:),'CO_{2}','{\it c} (mM)',x_DBL,b);
    nexttile; plotting2D(X,c_iHCO3(end,:),'HCO_{3}^{-}','{\it c} (mM)',x_DBL,b);
    nexttile; plotting2D(X,c_iCO3(end,:),'CO_{3}^{2-}','{\it c} (mM)',x_DBL,b);
    nexttile; plotting2D(X,c_iOH(end,:),'OH^{-}','{\it c} (mM)',x_DBL,b);
    nexttile; plotting2D(X,c_iH(end,:),'H^{+}','{\it c} (mM)',x_DBL,b);
    nexttile; plotting2D(X,pH(end,:),'pH','pH (-)',x_DBL,b);
    sgtitle('steady-state concentration profiles');
    fprintf('The profiles %.1f s after imposing %1.1f mA/cm2 at the cathode.\n', [t,j/10]);

    fprintf('\nConcentrations at the cathode surface after equilibration:\n');
    fprintf('\tCO2     %.5f mM\n',c_CO2(end,1));
    fprintf('\tHCO3(-) %.5f mM\n',c_iHCO3(end,1));
    fprintf('\tCO3(2-) %.5f mM\n',c_iCO3(end,1));
    fprintf('\tOH(-)   %.10f mM\n',c_iOH(end,1));
    fprintf('\tH(+)    %.10f mM\n',c_iH(end,1));
    fprintf('\tpH      %.2f   \n',pH(end,1));

    % Sanity check of confirming the pKw is adhered to in the equilibrium situation.
    figure(4); fs = 13;
    plot(X*1e6,-log10(c_iH(end,:).*c_iOH(end,:)/1e6),'color',b,LineWidth = 2); hold on;
    plot([0 x_DBL*1e6],[pKw(end,end) pKw(end,end)],'--k');
    set(gca,'Fontsize',fs);
    xlabel('{\it x} (μm)','FontSize',fs); xlim([0 x_DBL*1e6])
    ylabel('p{\it K}_{w} (-)', 'FontSize',fs); ylim([12 14.5])
    xtickformat('%.1f'); ytickformat('%.1f'); 
    title('water dissociation constant');
    fprintf('The calculated water dissociation constant compared to the theoretical value (dashed).\n');

    t4 = toc; fprintf('\nElapsed time for PART IV is %1.3f seconds.\n\n',t4);

    fprintf('\n***Total elapsed time for this script is %1.3f seconds.***',t1 + t2 + t3 + t4);
end

% Set up the governing eqns in the form accepted by the pdepe.m solver.
% The species are entered in the order CO2 (1), iHCO3 (2), iCO3 (3),
% iOH (4), and iH (5). If simple == 1, the simplified set of equations 
% is used.
function [c,f,s] = tdr(~,~,u,dudx,D,k,simple)
    c = [1; 1; 1; 1; 1];
    f = [D(1).CO2*dudx(1);
         D(1).iHCO3*dudx(2);
         D(1).iCO3*dudx(3);
         D(1).iOH*dudx(4);
         D(1).iH*dudx(5)];
    if simple == 1
        % ignores reaction terms based on order-of-magnitude analysis
        s = [-k(1).kb1b*u(1)*u(4) + k(1).kb1f*u(2) - k(1).ka1f*u(1) + k(1).ka1b*u(5)*u(2);
              k(1).kb1b*u(1)*u(4) - k(1).kb1f*u(2) - k(1).kb2b*u(2)*u(4) + k(1).kb2f*u(3) + k(1).ka1f*u(1) - k(1).ka1b*u(5)*u(2);
              k(1).kb2b*u(2)*u(4) - k(1).kb2f*u(3);
             -k(1).kb2b*u(2)*u(4) + k(1).kb2f*u(3);
             -k(1).ka2b*u(5)*u(3) + k(1).ka2f*u(2)];
    else
        % includes all reaction terms (acid/base components of (bi)carbonate equilibria and water dissociation)
        s = [-k(1).kb1b*u(1)*u(4) + k(1).kb1f*u(2) - k(1).ka1f*u(1)      + k(1).ka1b*u(5)*u(2);
              k(1).kb1b*u(1)*u(4) - k(1).kb1f*u(2) - k(1).kb2b*u(2)*u(4) + k(1).kb2f*u(3) + k(1).ka1f*u(1) - k(1).ka1b*u(5)*u(2) + k(1).ka2b*u(5)*u(3) - k(1).ka2f*u(2);
              k(1).kb2b*u(2)*u(4) - k(1).kb2f*u(3) - k(1).ka2b*u(5)*u(3) + k(1).ka2f*u(2);
             -k(1).kb1b*u(1)*u(4) + k(1).kb1f*u(2) - k(1).kb2b*u(2)*u(4) + k(1).kb2f*u(3) - k(1).kwb*u(4)*u(5) + k(1).kwf;
             -k(1).ka1b*u(5)*u(2) + k(1).ka1f*u(1) - k(1).ka2b*u(5)*u(3) + k(1).ka2f*u(2) - k(1).kwb*u(4)*u(5) + k(1).kwf];
    end
end

% Set up the TDR initial condition in the form accepted by the pdepe.m
% solver.
function u0 = tdr_ic(~,c_ic)
    u0 = [c_ic(1).CO2;
          c_ic(1).iHCO3;
          c_ic(1).iCO3;
          c_ic(1).iOH;
          c_ic(1).iH];
end

% Set up the TDR boundary conditions in the form accepted by the pdepe.m
% solver.
function [pl, ql, pr, qr] = tdr_bc(~,~,~,ur,~,c_ic,j,F,n)
    pl = [ j/F*(1/2*n(1) + 1/2*n(2) + 1/8*n(3) + 2/12*n(4) + 2/12*n(5) + 2/8*n(6));
           0;
           0;
          -j/F*(2/2*n(1) + 1/2*n(2) + 8/8*n(3) + 12/12*n(4) + 12/12*n(5) + 7/8*n(6) + 2/2*(1 - sum(n)));
           0];
    ql = [1; 1; 1; 1; 1];
    pr = [ur(1) - c_ic(1).CO2;
          ur(2) - c_ic(1).iHCO3;
          ur(3) - c_ic(1).iCO3;
          ur(4) - c_ic(1).iOH;
          ur(5) - c_ic(1).iH];
    qr = [0; 0; 0; 0; 0];
end

% Plotting the concentration profiles in space.
%  @X       vector of grid points in space
%  @c       vector of concentrations
%  @species name for the plot title
%  @yaxis   name for the y-axis
%  @x_DBL   right-hand limit for the x-axis
%  @colour  colour of the line
function plotting2D(X,c,species,yaxis,x_DBL,colour)
    fs = 13;
    plot(X*1e6,c,LineWidth = 2,Color = colour); hold on; title(species); set(gca,'Fontsize',fs);
    xlabel('{\it x} (μm)','FontSize',fs); ylabel(yaxis,'FontSize',fs);
    xtickformat('%.1f'); ytickformat('%.2f'); xlim([0 x_DBL*1e6]);
end

% Plotting the concentration profiles in time and space.
%  @X       vector of grid points in space
%  @T       vector of grid points in time
%  @c       matrix of concentrations
%  @species name for the plot title
%  @zaxis   name for the x-axis
function plotting3D(X,T,c,species,zaxis)
    fs = 13;
    surf(X.*1e6,T,c); title(species); set(gca,'Fontsize',fs);
    xlabel('{\it x} (μm)','FontSize',fs); ylabel('{\it t} (s)', 'FontSize',fs); zlabel(zaxis,'FontSize',fs);
    xtickformat('%.1f'); ytickformat('%.1f'); ztickformat('%.2f'); colorbar;
end