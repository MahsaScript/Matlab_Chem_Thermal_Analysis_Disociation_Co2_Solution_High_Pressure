% Equilibrium potentials of the HER and CO2RR calculated by
%
%    E_0_T = E_0 + dS_rxn/(nF)*(T - 298.15)
%
% The standard state equilibrium potentials are calculated as described
% in Appendix A and the entropy of reaction is calculated based on Hess'
% law using thermochemical data tables from the source on pages ...
%
% All equilibrium electrode potentials are on an SHE scale.

%  @T       scalar or vector of temperature in [K]
%  @returns structure as indicated below
function E_0 = nernst(T)
    % column 1: standard state equilibrium potential at T = 298.15 K in [V]
    % column 2: standard state entropy of reaction at T = 298.15 K in [J/mol/K]
    % column 3: electrons involved per extent of reaction in [-]
    % column 4:                equilibrium potential at T in [V]
    % column 5:                equilibrium potential at bulk conditions at T in [V]
    % column 6:                equilibrium potential at local conditions at T in [V]
    % field names preceded by 'i' indicate ionic species
    % note that columns 5 and 6 are evaluated during the simulation
    E_0 = struct('H2',     { 0.000, -31.1, 2,0,0,0},...
                 'CO',     {-0.104,-107.9, 2,0,0,0},...
                 'CH4',    { 0.169,0, 8,0,0,0},...
                 'iHCOO',  { 0.378,0, 2,0,0,0},... 
                 'C2H4',   { 0.079,0,12,0,0,0},...
                 'C2H5OH', { 0.084,0,12,0,0,0},...
                 'iCH3COO',{ 0.150,0, 8,0,0,0});     

    % overwrite the previous structure with corrected OCP values
    fields = fieldnames(E_0);
    for i  = 1:numel(fields)
        E_0(4).(fields{i}) = E_0(1).(fields{i}) + E_0(2).(fields{i})/(E_0(3).(fields{i})*96485).*(T' - 298.15);
    end
end