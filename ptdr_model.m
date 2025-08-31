% This program may be copied and used in unaltered form. Changes may be
% made to the contents at any point for personal use. Redistribution of
% this program in both altered and unaltered form is subject to the
% conditions that it remains free of charge and that this notice is
% included upfront. No warranties are given, expressed or implied, that
% this program is free of errors and bugs. Neither do we guarantee that
% this program meets the requirements of your particular application. The
% author disclaims all liability for direct or consequential damages as a
% result from use of this program.

function ptdr_model(T,P,j,n,electrolyte,c,e,l,k_m,L,n_xCL,n_xDBL,t,dt)


% This model uses the pdepe.m solver to find numerical solutions to the
% governing transport equations of the form
%
%   d/dt(c) = D*d^2/dx^2(c) + R_V 
%
% which describes the transient reaction-diffusion equation for species c
% with volumetric reaction terms R_V (which couples the various equations).
% This is a 1D system describing the porous catalyst layer between x = 0
% and x = l and diffusion boundary layer between x = l and x = L.
%
% Note that the concentrations returned by the solver are in units of
% [mol/m3] (i.e. [mM]) to keep in line with the dimensionality of the
% governing equation. For equilibrium and rate coefficients, their values
% are outputted in the console window in dimensions involving [mol/L] for
% ease of comparison with literature values, after which they are converted
% to [mol/m3]-equivalent values and used as such in the solution procedure.

    clc; close all;
    format compact; format short;
    fprintf('***1D Transient Diffusion-Reaction Model for a Porous Electrode***\n\n')

    % cast input parameters in other units
    j             = j*10;                      % current density [A/m2]
    T             = 273.15 + T;                % temperature [K]
    x_CL          = l/1e6;                     % catalyst layer thickness [m]
    x_DBL         = L/1e6;                     % diffusion layer thickness [m]
    % values of other natural constants
    F             = 96485;                     % Faraday's constant [C/mol]
    % plot colours
    b = [0 0.4470 0.7410];
    o = [0.8500 0.3250 0.0980];
    y = [0.9290 0.6940 0.1250];
    g = [0.4660 0.6740 0.1880];

    tic
    % (1) Initialize parameters relevant to the model.
    fprintf('\n---PART I: Thermodynamic and Kinetic Parameters---\n\n');
    % viscosity of KHCO3 electrolyte in [Pa s] as function of molality
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

    D_CL = diffusivity(T,u,0,e);
    fprintf('Diffusion coefficients of species in the catalyst layer (porosity- and tortuosity-corrected):\n')
    fprintf('\t D(CO2)     = %2.3e m2/s\n',D_CL(1).CO2);
    fprintf('\t D(NO3(-))  = %2.3e m2/s\n',D_CL(1).iNO3);
    fprintf('\t D(HCO3(-)) = %2.3e m2/s\n',D_CL(1).iHCO3);
    fprintf('\t D(CO3(2-)) = %2.3e m2/s\n',D_CL(1).iCO3);        
    fprintf('\t D(OH(-))   = %2.3e m2/s\n',D_CL(1).iOH);
    fprintf('\t D(H(+))    = %2.3e m2/s\n',D_CL(1).iH);
    fprintf('\t D(NH4(+))  = %2.3e m2/s\n',D_CL(1).iNH4);
    D_DBL = diffusivity(T,u,0,1);
    fprintf('Diffusion coefficients of species in the diffusion layer:\n')
    fprintf('\t D(CO2)     = %2.3e m2/s\n',D_DBL(1).CO2);
    fprintf('\t D(NO3(-))  = %2.3e m2/s\n',D_DBL(1).iNO3);
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
    fprintf('Corrected for temperature and viscosity.\n\n');

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
    fprintf('\t kw>         = %2.3e M/s\t\t kw<        = %2.3e /M/s\n',[k.kwf, k.kwb]);
    fprintf('Corrected for temperature and viscosity.\n\n');
    % rewrite all kinetic coefficients in units involving [m3]
    k(1).ka1b = k(1).ka1b/1e3; k(1).ka2b = k(1).ka2b/1e3; k(1).kb1b = k(1).kb1b/1e3;
    k(1).kb2b = k(1).kb2b/1e3; k(1).kwf  = k(1).kwf*1e3;  k(1).kwb  = k(1).kwb/1e3;
    t1 = toc; fprintf('\nElapsed time for PART I is %1.3f seconds.\n\n',t1);

    tic
    % (2) Calculation of concentrations in the bulk solution.
    fprintf('\n---PART II: Bulk Equilibrium Concentrations---\n');
    [~,~,~,~] = bjerrum(Kw,Kb1,Kb2);
    title('Bjerrum plot of the carbonate system');
    fprintf('\nThe fraction of CO2/HCO3(-)/CO3(2-) of the Dissolved Inorganic Content (DIC) as a function of pH.\n');
    if electrolyte == 'KOH'
        % using a root finding algorithm to find the initial conditions
        p_coeff   = [2*3e-3*c/1e3 (1 + 3e-3*c/1e3/Kb1)*Kb1*Kb2 -c/1e3*Kb1*Kb2 -Kw*Kb1*Kb2];
        r         = roots(p_coeff);
        c_ic      = struct('CO2',   {0},...
                           'iHCO3', {0},...
                           'iCO3',  {0},...
                           'iOH',   {c},...
                           'iH',    {Kw*1e6/c});
        pH_i      = pKw + log10(c_ic(1).iOH/1e3);
    elseif electrolyte == 'KHCO3'
        % using a root finding algorithm to find the initial conditions
        p_coeff   = [2*3e-3*c/1e3 (1 + 3e-3*c/1e3/Kb1)*Kb1*Kb2 -c/1e3*Kb1*Kb2 -Kw*Kb1*Kb2];
        r         = roots(p_coeff);
        c_ic      = struct('CO2',   {3e-3*c},...
                           'iHCO3', {3e-3*c*r(2)/Kb1},...
                           'iCO3',  {3e-3*c*r(2)/Kb1*r(2)/Kb2},...
                           'iOH',   {r(2)*1e3},...
                           'iH',    {Kw/r(2)*1e3});
        pH_i      = -log10(c_ic.iH/1e3);
    elseif electrolyte == 'KCl'
        pH_i      = pKw/2;
        c_ic      = struct('CO2',   {0},...
                           'iHCO3', {0},...
                           'iCO3',  {0},...
                           'iOH',   {1e3*sqrt(Kw)},...
                           'iH',    {1e3*sqrt(Kw)});
    elseif electrolyte == 'K2CO3'
        pH_i      = pKw + log10(sqrt(c/1e3*Kb2));
        c_ic      = struct('CO2',   {0},...
                           'iHCO3', {1e3*sqrt(c/1e3*Kb2)},...
                           'iCO3',  {c},...
                           'iOH',   {1e3*10^-(pKw - pH_i)},...
                           'iH',    {1e3*10^(-pH_i)});
    end
    fprintf('\nThe initial conditions of the electrolyte are:\n');
    fprintf('\tCO2     %.5f mM\n',c_ic(1).CO2);
    fprintf('\tOH(-)   %.5f mM\n',c_ic(1).iOH);
    fprintf('\tHCO3(-) %.5f mM\n',c_ic(1).iHCO3);
    fprintf('\tCO3(2-) %.5f mM\n',c_ic(1).iCO3);
    fprintf('\tH(+)    %.5f mM\n',c_ic(1).iH);
    fprintf('\tpH      %.2f \n',pH_i);
    t2 = toc; fprintf('\nElapsed time for PART II is %1.3f seconds.\n\n',t2);

    tic
    % (3) Setting up the system of matrix-vector equations for the model.
    fprintf('\n---PART III: Concentration Profiles---\n');
    fprintf('\nModel parameters:\n');
    fprintf('\tcatalyst layer thickness               %.2f um\n',x_CL*1e6);
    fprintf('\tdiffusion boundary layer thickness     %.2f um\n',x_DBL*1e6);
    fprintf('\tlogarithmic grid spacing in CL         %.0f points\n',n_xCL);
    fprintf('\tlogarithmic grid spacing in DBL        %.0f points\n',n_xDBL);
    fprintf('\ttotal modelling time                   %.5f s\n',t);
    fprintf('\ttime step size                         %.5f s\n',dt);

    % discretize the catalyst layer domain logarithmically
    dx_CL        = x_CL/(n_xCL - 1);
    x_CL_l       = horzcat(0,logspace(log10(dx_CL),log10(x_CL/2),(n_xCL - 1)/2));
    x_CL_r       = x_CL - x_CL_l; x_CL_r = x_CL_r(1:end - 1); x_CL_r = flip(x_CL_r);
    X_CL         = horzcat(x_CL_l,x_CL_r);
    % discretize the diffusion layer domain logarithmically
    dx_DBL       = 1e-1*x_DBL/(n_xDBL - 1);
    x_DBL_l      = horzcat(0,logspace(log10(dx_DBL),log10(x_DBL/2),(n_xDBL - 1)/2));
    x_DBL_r      = x_DBL - x_DBL_l; x_DBL_r = x_DBL_r(1:end - 1); x_DBL_r = flip(x_DBL_r);
    X_DBL        = x_CL + horzcat(x_DBL_l,x_DBL_r);
    X            = horzcat(X_CL,X_DBL(2:end));
    % discretize the time domain linearly
    N            = [0:dt:t];

    cCO2 = zeros(4,n_xCL + n_xDBL - 1);
    cpH  = zeros(4,n_xCL + n_xDBL - 1);
    tCE = zeros(2,length(j));
    for i = 1:length(j)
        % The volumetric reaction rate for CO2 consumption and OH(-) generation.
        R_CO2 =  j(i)/(F*x_CL*e)*(1/2*n(i,1) + 1/2*n(i,2) + 1/8*n(i,3) + 2/12*n(i,4) + 2/12*n(i,5) + 2/8*n(i,6));
        R_iOH = -j(i)/(F*x_CL*e)*(2/2*n(i,1) + 1/2*n(i,2) + 8/8*n(i,3) + 12/12*n(i,4) + 12/12*n(i,5) + 7/8*n(i,6) + 2/2*(1 - sum(n(i,:))));
        fprintf('The averaged volumetric reaction rates in the catalyst layer evaluate to:\n');
        fprintf('\tR_CO2     %.5f mol/m3/s\n',R_CO2);
        fprintf('\tR_OH(-)   %.5f mol/m3/s\n',R_iOH);
        fprintf('Note that this value is normalized to the pore volume of the catalyst layer!\n');
        
        % solve the system of equations
        param_tdr    = @(x,t,u,dudx)ptdr(x,t,u,dudx,D_CL,D_DBL,k,R_CO2,R_iOH,x_CL,electrolyte);
        param_tdr_ic = @(x)tdr_ic(x,c_ic,c_sat,x_CL,electrolyte);
        param_tdr_bc = @(xl,ul,xr,ur,t)tdr_bc(xl,ul,xr,ur,t,c_ic,c_sat,k_m,electrolyte);
        sol     = pdepe(0,param_tdr,param_tdr_ic,param_tdr_bc,X,N);
        % dim 1: time T; dim 2: position X; dim 3: species type
        c_CO2   = sol(:,:,1);
        if strcmp(electrolyte,'KHCO3')
            c_iHCO3 = exp(sol(:,:,2));
            c_iCO3  = exp(sol(:,:,3));
        elseif strcmp(electrolyte,'KOH')
            c_iHCO3 = sol(:,:,2);
            c_iCO3  = sol(:,:,3);
        end
        c_iOH   = exp(sol(:,:,4)); 
        c_iH    = exp(sol(:,:,5));
        pKw     = -log10(c_iH.*c_iOH/1e6);
        pH      = pKw + log10(c_iOH/1e3);
    
        % checking if the CO2 concentration is negative.
        if c_CO2(end,n_xCL) < 0
            fprintf('\nWarning: The partial current density exceeds the maximum supportable for CO2 reduction for the combination n_CO = %.2f and j = %.1f mA/cm2.\n',[n(1),j(i)/10])
        end

        fprintf('\nConcentrations at the GDL-CL interface after %.3f s:\n', [t]);
        fprintf('\tCO2     %.5f mM\n',c_CO2(end,1));
        fprintf('\tHCO3(-) %.5f mM\n',c_iHCO3(end,1));
        fprintf('\tCO3(2-) %.5f mM\n',c_iCO3(end,1));
        fprintf('\tOH(-)   %.5f mM\n',c_iOH(end,1));
        fprintf('\tH(+)    %.5f mM\n',c_iH(end,1));
        fprintf('\tpH      %.2f   \n',pH(end,1));
        fprintf('\nConcentrations at the CL-DBL interface after %.3f s:\n', [t]);
        fprintf('\tCO2     %.5f mM\n',c_CO2(end,n_xCL));
        fprintf('\tHCO3(-) %.5f mM\n',c_iHCO3(end,n_xCL));
        fprintf('\tCO3(2-) %.5f mM\n',c_iCO3(end,n_xCL));
        fprintf('\tOH(-)   %.5f mM\n',c_iOH(end,n_xCL));
        fprintf('\tH(+)    %.5f mM\n',c_iH(end,n_xCL));
        fprintf('\tpH      %.2f   \n',pH(end,n_xCL));
        CE = carbonEff(c_CO2,D_CL,X,R_CO2,x_CL,n_xCL);
        fprintf('\nThe carbon efficiencies according to the below definitions:\n')
        fprintf('\tcarbon efficiency #1      = %.3f %.', [CE(1)]);
        fprintf('\n\tcarbon efficiency #2      = %.3f %.', [CE(2)]);
        fprintf('\n\toverall carbon efficiency = %.3f %.', [CE(1).*CE(2)]);
        fprintf('\nDef. #1 = fraction of CO2 consumed in the CL compared to that entering the system at x = 0');
        fprintf('\nDef. #2 = fraction of CO2 utilized for CO2RR in the CL compared to all consumptions in the CL');

        cCO2(i,:) = c_CO2(end,:);
        cpH(i,:)  = pH(end,:);
        tCE(1,i) = CE(1);
        tCE(2,i) = CE(2);
    end
    t3 = toc; fprintf('\nElapsed time for PART III is %1.3f seconds.\n\n',t3);

    % (4) Preparing plots of the model solutions.
    tic

    figure(2);
    tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
    nexttile; plotting3D(X,N,c_CO2,'CO_{2}');
    nexttile; plotting3D(X,N,c_iHCO3,'HCO_{3}^{-}');
    nexttile; plotting3D(X,N,c_iCO3,'CO_{3}^{2-}');
    nexttile; plotting3D(X,N,c_iOH,'OH^{-}');
    nexttile; plotting3D(X,N,c_iH,'H^{+}');
    nexttile; plotting3D(X,N,pH,'pH');
    sgtitle('transient concentration profiles');
    fprintf('The transients upon imposing %1.1f mA/cm2 at the cathode.\n', j/10);

    figure(3);
    b = [0.8500 0.3250 0.0980];
    tiledlayout(3,2,'TileSpacing','Compact','Padding','Compact');
    nexttile; plotting2D(X,c_CO2(end,:),'CO_{2}',x_CL,x_DBL,b);
    nexttile; plotting2D(X,c_iHCO3(end,:),'HCO_{3}^{-}',x_CL,x_DBL,b);
    nexttile; plotting2D(X,c_iCO3(end,:),'CO_{3}^{2-}',x_CL,x_DBL,b);
    nexttile; plotting2D(X,c_iOH(end,:),'OH^{-}',x_CL,x_DBL,b);
    nexttile; plotting2D(X,c_iH(end,:),'H^{+}',x_CL,x_DBL,b);
    nexttile; plotting2D(X,pH(end,:),'pH',x_CL,x_DBL,b);
    sgtitle('steady-state concentration profiles');
    fprintf('The profiles %.1f s after imposing %1.1f mA/cm2 at the cathode.\n', [t,j/10]);

    figure(4); fs = 13;
    plot(X*1e6,-log10(c_iH(end,:).*c_iOH(end,:)/1e6),'color',b,LineWidth = 2); hold on;
    plot([0 x_DBL*1e6],[pKw(end,end) pKw(end,end)],'--k');
    set(gca,'Fontsize',fs);
    xlabel('{\it x} (μm)','FontSize',fs); xlim([0 (x_DBL+x_CL)*1e6])
    ylabel('p{\it K}_{w} (-)', 'FontSize',fs); ylim([12 15])
    xtickformat('%.1f'); ytickformat('%.1f'); 
    title('water dissociation constant');
    fprintf('The calculated water dissociation constant compared to the theoretical value (dashed).\n');

    t4 = toc; fprintf('\nElapsed time for PART IV is %1.3f seconds.\n\n',t4);

%     figure(6);
%     % plot colours
%     b = [0 0.4470 0.7410];
%     o = [0.8500 0.3250 0.0980];
%     y = [0.9290 0.6940 0.1250];
%     g = [0.4660 0.6740 0.1880];
%     fs = 13;
%     plot(X*1e6,cCO2(1,:),LineWidth = 2,Color = b); hold on;
%     plot(X*1e6,cCO2(2,:),LineWidth = 2,Color = o); hold on;
%     plot(X*1e6,cCO2(3,:),LineWidth = 2,Color = y); hold on;
%     plot(X*1e6,cCO2(4,:),LineWidth = 2,Color = g);
%     set(gca,'Fontsize',fs);
%     xlabel('{\it x} (μm)','FontSize',fs); ylabel('{\it c} (mM)','FontSize',fs);
%     xtickformat('%.0f'); ytickformat('%.0f'); xlim([0 (x_CL+x_DBL)*1e6]);
% 
%     figure(7);
%     plot(X*1e6,cpH(1,:),LineWidth = 2,Color = b); hold on;
%     plot(X*1e6,cpH(2,:),LineWidth = 2,Color = o); hold on;
%     plot(X*1e6,cpH(3,:),LineWidth = 2,Color = y); hold on;
%     plot(X*1e6,cpH(4,:),LineWidth = 2,Color = g);
%     set(gca,'Fontsize',fs);
%     xlabel('{\it x} (μm)','FontSize',fs); ylabel('pH (-)','FontSize',fs);
%     xtickformat('%.0f'); ytickformat('%.0f'); xlim([0 (x_CL+x_DBL)*1e6]);
%     legend('-10 mA cm^{-2}','-50 mA cm^{-2}','-100 mA cm^{-2}','-200 mA cm^{-2}')

    fprintf('\n***Total elapsed time for this script is %1.3f seconds.***',t1 + t2 + t3 + t4);
end

% Set up the governing eqns in the form accepted by the pdepe.m solver.
% The species are entered in the order CO2 (1), iHCO3 (2), iCO3 (3),
% iOH (4), and iH (5).
function [c,f,s] = ptdr(x,~,u,dudx,D_CL,D_DBL,k,R_CO2,R_iOH,x_CL,electrolyte)
    c = [1; 1; 1; 1; 1];
    if strcmp(electrolyte,'KHCO3')
        if x <= x_CL
            f = [D_CL(1).CO2*dudx(1);
                 exp(u(2))*D_CL(1).iHCO3*dudx(2);
                 exp(u(3))*D_CL(1).iCO3*dudx(3);
                 exp(u(4))*D_CL(1).iOH*dudx(4);
                 exp(u(5))*D_CL(1).iH*dudx(5)];
            s = [-k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*exp(u(2)) - k(1).ka1f*u(1) + k(1).ka1b*exp(u(2))*exp(u(5)) + R_CO2;
                  k(1).kb1b*u(1)*exp(u(4)) - k(1).kb1f*exp(u(2)) - k(1).kb2b*exp(u(2))*exp(u(4)) + k(1).kb2f*exp(u(3)) + k(1).ka1f*u(1) - k(1).ka1b*exp(u(2))*exp(u(5)) + k(1).ka2b*exp(u(3))*exp(u(5)) - k(1).ka2f*exp(u(2));
                  k(1).kb2b*exp(u(2))*exp(u(4)) - k(1).kb2f*exp(u(3)) - k(1).ka2b*exp(u(3))*exp(u(5)) + k(1).ka2f*exp(u(2));
                 -k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*exp(u(2)) - k(1).kb2b*exp(u(2))*exp(u(4)) + k(1).kb2f*exp(u(3)) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf + R_iOH;
                 -k(1).ka1b*exp(u(2))*exp(u(5)) + k(1).ka1f*u(1) - k(1).ka2b*exp(u(3))*exp(u(5)) + k(1).ka2f*exp(u(2)) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf];
        else
            f = [D_DBL(1).CO2*dudx(1);
                 exp(u(2))*D_DBL(1).iHCO3*dudx(2);
                 exp(u(3))*D_DBL(1).iCO3*dudx(3);
                 exp(u(4))*D_DBL(1).iOH*dudx(4);
                 exp(u(5))*D_DBL(1).iH*dudx(5)];
            s = [-k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*exp(u(2)) - k(1).ka1f*u(1) + k(1).ka1b*exp(u(2))*exp(u(5));
                  k(1).kb1b*u(1)*exp(u(4)) - k(1).kb1f*exp(u(2)) - k(1).kb2b*exp(u(2))*exp(u(4)) + k(1).kb2f*exp(u(3)) + k(1).ka1f*u(1) - k(1).ka1b*exp(u(2))*exp(u(5)) + k(1).ka2b*exp(u(3))*exp(u(5)) - k(1).ka2f*exp(u(2));
                  k(1).kb2b*exp(u(2))*exp(u(4)) - k(1).kb2f*exp(u(3)) - k(1).ka2b*exp(u(3))*exp(u(5)) + k(1).ka2f*exp(u(2));
                 -k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*exp(u(2)) - k(1).kb2b*exp(u(2))*exp(u(4)) + k(1).kb2f*exp(u(3)) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf;
                 -k(1).ka1b*exp(u(2))*exp(u(5)) + k(1).ka1f*u(1) - k(1).ka2b*exp(u(3))*exp(u(5)) + k(1).ka2f*exp(u(2)) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf];
        end
    elseif strcmp(electrolyte,'KOH')
        if x <= x_CL
            f = [D_CL(1).CO2*dudx(1);
                 D_CL(1).iHCO3*dudx(2);
                 D_CL(1).iCO3*dudx(3);
                 exp(u(4))*D_CL(1).iOH*dudx(4);
                 exp(u(5))*D_CL(1).iH*dudx(5)];
            s = [-k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*u(2) - k(1).ka1f*u(1)           + k(1).ka1b*exp(u(5))*u(2) + R_CO2;
                  k(1).kb1b*u(1)*exp(u(4)) - k(1).kb1f*u(2) - k(1).kb2b*u(2)*exp(u(4)) + k(1).kb2f*u(3) + k(1).ka1f*u(1) - k(1).ka1b*exp(u(5))*u(2) + k(1).ka2b*exp(u(5))*u(3) - k(1).ka2f*u(2);
                  k(1).kb2b*u(2)*exp(u(4)) - k(1).kb2f*u(3) - k(1).ka2b*exp(u(5))*u(3) + k(1).ka2f*u(2);
                 -k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*u(2) - k(1).kb2b*u(2)*exp(u(4)) + k(1).kb2f*u(3) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf + R_iOH;
                 -k(1).ka1b*exp(u(5))*u(2) + k(1).ka1f*u(1) - k(1).ka2b*exp(u(5))*u(3) + k(1).ka2f*u(2) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf];
        else
            f = [D_DBL(1).CO2*dudx(1);
                 D_DBL(1).iHCO3*dudx(2);
                 D_DBL(1).iCO3*dudx(3);
                 exp(u(4))*D_DBL(1).iOH*dudx(4);
                 exp(u(5))*D_DBL(1).iH*dudx(5)];
            s = [-k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*u(2) - k(1).ka1f*u(1)           + k(1).ka1b*exp(u(5))*u(2);
                  k(1).kb1b*u(1)*exp(u(4)) - k(1).kb1f*u(2) - k(1).kb2b*u(2)*exp(u(4)) + k(1).kb2f*u(3) + k(1).ka1f*u(1) - k(1).ka1b*exp(u(5))*u(2) + k(1).ka2b*exp(u(5))*u(3) - k(1).ka2f*u(2);
                  k(1).kb2b*u(2)*exp(u(4)) - k(1).kb2f*u(3) - k(1).ka2b*exp(u(5))*u(3) + k(1).ka2f*u(2);
                 -k(1).kb1b*u(1)*exp(u(4)) + k(1).kb1f*u(2) - k(1).kb2b*u(2)*exp(u(4)) + k(1).kb2f*u(3) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf;
                 -k(1).ka1b*exp(u(5))*u(2) + k(1).ka1f*u(1) - k(1).ka2b*exp(u(5))*u(3) + k(1).ka2f*u(2) - k(1).kwb*exp(u(4))*exp(u(5)) + k(1).kwf];
        end
    end
end

% Set up the TDR initial condition in the form accepted by the pdepe.m
% solver.
function u0 = tdr_ic(x,c_ic,c_sat,l_CL,electrolyte)
    if strcmp(electrolyte,'KHCO3')
        if x <= l_CL
            u0 = [c_sat;
                  log(c_ic(1).iHCO3);
                  log(c_ic(1).iCO3);
                  log(c_ic(1).iOH);
                  log(c_ic(1).iH)];
        else
            u0 = [c_ic(1).CO2;
                  log(c_ic(1).iHCO3);
                  log(c_ic(1).iCO3);
                  log(c_ic(1).iOH);
                  log(c_ic(1).iH)];
        end
    elseif strcmp(electrolyte,'KOH')
        if x <= l_CL
            u0 = [c_sat;
                  c_ic(1).iHCO3;
                  c_ic(1).iCO3;
                  log(c_ic(1).iOH);
                  log(c_ic(1).iH)];
        else
            u0 = [c_ic(1).CO2;
                  c_ic(1).iHCO3;
                  c_ic(1).iCO3;
                  log(c_ic(1).iOH);
                  log(c_ic(1).iH)];
        end
    end
end

% Set up the TDR boundary conditions in the form accepted by the pdepe.m
% solver.
function [pl, ql, pr, qr] = tdr_bc(~,ul,~,ur,~,c_ic,c_sat,k_m,electrolyte)
    pl = [k_m*(c_sat - ul(1));
          0;
          0;
          0;
          0];
    ql = [1; 1; 1; 1; 1];
    if strcmp(electrolyte,'KHCO3')
        pr = [ur(1) - c_ic(1).CO2;
              ur(2) - log(c_ic(1).iHCO3);
              ur(3) - log(c_ic(1).iCO3);
              ur(4) - log(c_ic(1).iOH);
              ur(5) - log(c_ic(1).iH)];
    elseif strcmp(electrolyte,'KOH')
        pr = [ur(1) - c_ic(1).CO2;
              ur(2) - c_ic(1).iHCO3;
              ur(3) - c_ic(1).iCO3;
              ur(4) - log(c_ic(1).iOH);
              ur(5) - log(c_ic(1).iH)];
    end
    qr = [0; 0; 0; 0; 0];
end

% Plotting the concentration profiles in space.
%  @X       vector of grid points in space
%  @c       vector of concentrations
%  @species name for the plot title
%  @x_L     right-hand limit for the x-axis
%  @colour  colour of the line
function plotting2D(X,c,species,x_CL,x_DBL,colour)
    fs = 13;
    plot(X*1e6,c,LineWidth = 2,Color = colour); hold on; 
    plot([x_CL,x_CL]*1e6,[0,1.2*max(c)],LineWidth = 2,LineStyle = ':',Color = 'k'); title(species); set(gca,'Fontsize',fs);
    xlabel('{\it x} (μm)','FontSize',fs); ylabel('{\it c} (mM)','FontSize',fs);
    xtickformat('%.1f'); ytickformat('%.2f'); xlim([0 (x_CL+x_DBL)*1e6]);
end

% Plotting the concentration profiles in time and space.
%  @X       vector of grid points in space
%  @T       vector of grid points in time
%  @c       matrix of concentrations
%  @species name for the plot title
function plotting3D(X,T,c,species)
    fs = 13;
    surf(X.*1e6,T,c); title(species); set(gca,'Fontsize',fs);
    xlabel('{\it x} (μm)','FontSize',fs); ylabel('{\it t} (s)', 'FontSize',fs); zlabel('{\it c} (mM)','FontSize',fs);
    xtickformat('%.1f'); ytickformat('%.1f'); ztickformat('%.2f'); colorbar;
end

% Calculating carbon efficiency (CE) by two definitions: (1) the fraction 
% of CO2 escaping into the electrolyte versus that entering the CL; (2) the
% fraction of CO2 being used in electrochemical reactions versus that being
% lost to (bi)carbonate equilibria.
%  @c_CO2   concentration matrix of CO2 from the pdepe solver
%  @D_CL    diffusion coefficient of CO2 in the CL in [m2/s]
%  @X       vector of node locations of discretized domain in [m]
%  @R_CO2   reaction rate of CO2 averaged over the entire CL volume in [mol/m3/s]
%  @x_CL    thickness of the CL in [m]
%  @nx_CL   number of grid points in the CL in [-]
%  @returns vector of DBL and CL efficiency resp. in [-]
function CE = carbonEff(c_CO2,D_CL,X,R_CO2,x_CL,n_xCL)
    Jin   = -D_CL(1).CO2*(c_CO2(end,2) - c_CO2(end,1))/(X(2)-X(1));                           % [mol/m2/s] 
    Jout  = -D_CL(1).CO2*(c_CO2(end,n_xCL + 1) - c_CO2(end,n_xCL))/(X(n_xCL + 1) - X(n_xCL)); % [mol/m2/s]
    R     = R_CO2*x_CL;                                                                       % [mol/m2/s]
    CE    = [1 - Jout/Jin; -R/(Jin - Jout)];
end