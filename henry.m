% Calculates the Henry solubility of the given species based on pressure
% (bar) and the Henry's Law solubility constant H [mol/m3/Pa] defined at
% 298.15 K (see section 2.4.1 in the source).

%  @T       scalar or vector of temperature in [K]
%  @P       scalar or vector of pressure in [bar]
%  @species string element (i.e. 'CO2', 'H2', etc...)
%  @returns scalar, vector, or matrix of CO2-saturated concentration in [mol/L]
function c_sat = henry(T,P,species)

    % Model parameters to describe temperature dependence of the Henry con-
    % stant. The first parameter is H_cp [mol/m3/Pa] and the second parame-
    % ter is d/d(1/T)[(ln(H_cp)] [K], which is equivalent to -ΔH_sol/R and
    % is used according to Eqn. (19) in the source:
    %
    %    H(T) = H_0*exp(-ΔH_sol/R*(1/T - 1/298.15))
    %
    H     = struct('H2',    { 7.8e-6,  530},...
                   'N2',    { 6.4e-6, 1600},...
                   'O2',    { 1.2e-5, 1700},...
                   'CO',    { 9.7e-6, 1300},... 
                   'CO2',   { 3.3e-4, 2400},... 
                   'CH4',   { 1.4e-5, 1900},...
                   'C2H4',  { 5.9e-5, 2200},...
                   'C2H6',  { 1.9e-5, 2400},...
                   'C3H6',  { 5.5e-5, 2800},...
                   'C3H8',  { 1.5e-5, 2700});
    P     = 1e5*P';
    H_cp  = H(1).(species)*exp(H(2).(species)*(1./T - 1/298.15));
    % in case vectors of temperature and pressure were specified
    %    T-dependence across a row of the matrix
    %    P-dependence across a column of the matrix
    c_sat = P*H_cp/1e3;
end