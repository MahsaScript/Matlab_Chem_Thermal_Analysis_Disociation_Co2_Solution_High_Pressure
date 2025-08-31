% Current densities according to Butler-Volmer kinetics including kinetic
% and transport polarization. The species currently supported by this are 
% H2, CO, CH3COO(-), CH4, C2H4, C2H5OH, and CH3COO(-). Charge transfer
% coefficient and exchange current density are to be specified by the user
% as vector covering all products in the order just mentioned. 
%
%  @j_0     vector of exchange current density in [A/m2]
%  @a       vector of charge transfer coefficient in [-]
%  @T       temperature in [K]
%  @E       applied potential on SHE scale in [V]
%  @E_0     equilibrium potential on SHE scale in [V]
%  @c       surface concentration in [mol/m3]
%  @c_0     equilibrium concentration in [mol/m3]
%  @returns structure with current density of product in [A/m2]
%                          total polarization in [V]
%                          transport polarization in [V]
%                          kinetic polarization in [V]
function kin = butlerVolmer(j_0,a,T,E,E_0,c,c_0)
    R = 8.3145; % [J/mol/K]
    F = 96485;  % [C/mol]
    % column 1: partial current density in [A/m2]
    % column 2: total polarization in [V]
    % column 3: transport polarization (T) in [V]
    % column 4: transport polarization (K) in [V]
    % column 5: kinetic polarization in [V]
    % field names preceded by 'i' indicate ionic species
    kin       = struct('H2',     {0,0,0,0},...
                       'CO',     {0,0,0,0},...
                       'CH4',    {0,0,0,0},...
                       'iHCOO',  {0,0,0,0},... 
                       'C2H4',   {0,0,0,0},...
                       'C2H5OH', {0,0,0,0},...
                       'iCH3COO',{0,0,0,0});  

    kin(1).H2    = -j_0(1)*exp(-2*a(1)*F/(R*T)*(E - E_0(6).H2));
    kin(2).H2    = E - E_0(5).H2;
    kin(3).H2    = E_0(6).H2 - E_0(5).H2;
    kin(4).H2    = 0;
    kin(5).H2    = R*T/(2*a(2)*F)*log(j_0(2)./kin(1).H2);
    kin(1).CO    = -j_0(2)*c(2)/c_0(2)*exp(-2*a(2)*F/(R*T)*(E - E_0(6).CO));
    kin(2).CO    = E - E_0(5).CO;
    kin(3).CO    = E_0(6).CO - E_0(5).CO;
    kin(4).CO    = R*T/(2*a(2)*F)*log(c(2)/c_0(2));
    kin(5).CO    = R*T/(2*a(2)*F)*log(j_0(2)./kin(1).CO);
    kin(1).iHCOO = -j_0(3)*c(2)/c_0(2)*c_0(1)/c(1)*exp(-2*a(3)*F/(R*T)*(E - E_0(6).iHCOO));
    kin(2).iHCOO = E - E_0(5).iHCOO;
    kin(3).iHCOO = E_0(6).iHCOO - E_0(5).iHCOO;
    kin(4).iHCOO = R*T/(2*a(3)*F)*log(c(2)/c_0(2)*c_0(1)/c(1));
    kin(5).iHCOO = R*T/(2*a(3)*F)*log(j_0(3)./kin(1).iHCOO);
    kin(1).CH4   = -j_0(4)*c(2)/c_0(2)*(c_0(1)/c(1))^2*exp(-8*a(4)*F/(R*T)*(E - E_0(6).CH4));
    kin(2).CH4   = E - E_0(5).CH4;
    kin(3).CH4   = E_0(6).CH4 - E_0(5).CH4;
    kin(4).CH4   = R*T/(8*a(4)*F)*log(c(2)/c_0(2)*(c_0(1)/c(1))^2);
    kin(5).CH4   = R*T/(8*a(4)*F)*log(j_0(4)./kin(1).CH4);
    kin(1).C2H4  = -j_0(5)*(c(2)/c_0(2))^2*exp(-12*a(5)*F/(R*T)*(E - E_0(6).C2H4));
    kin(2).C2H4  = E - E_0(5).C2H4;
    kin(3).C2H4  = E_0(6).C2H4 - E_0(5).C2H4;
    kin(4).C2H4  = R*T/(6*a(5)*F)*log(c(2)/c_0(2));
    kin(5).C2H4  = R*T/(6*a(5)*F)*log(j_0(5)./kin(1).C2H4);
end