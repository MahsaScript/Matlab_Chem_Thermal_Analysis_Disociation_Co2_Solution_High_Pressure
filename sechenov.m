% Calculation of the temperature and/or ion dependent solubility of gaseous
% species according to Eqn. (1) from the source:
%
%    log(c_G,0/c_G) = K*c_s
%

%  @c       scalar or vector of salt/ion concentration in [mol/L]
%  @T       scalar or vector of temperature in [K]
%  @P       scalar pressure in [bar]
%  @species string element (i.e. 'CO2', 'H2', etc...)
%  @salt    string element (i.e. 'KHCO3', 'KOH', etc...)
%  @returns scalar, vector, or matrix of CO2-saturated concentration in [mol/L]
function c_sat_S = sechenov(c,T,P,species,salt)

    % Model parameters to describe temperature dependence in the Sechenov
    % relation. The first parameter is h_G,0 [L/mol] and the second para-
    % meter is h_T [L/mol/K]. The gas-specific model parameter is calcul-
    % ated according to Eqn. (4) in the source:
    %
    %    h_G = h_G,0 + h_T(T - 298.15)
    %
    h   = struct('H2',    {-0.0218, -0.299/1e3},... % 273-353 K
                 'N2',    {-0.0010, -0.605/1e3},... % 278-345 K
                 'O2',    { 0.0000, -0.334/1e3},... % 273-353 K
                 'CO',    { 0.0000,  0.000    },... % data not available
                 'CO2',   {-0.0172, -0.338/1e3},... % 273-313 K
                 'CH4',   { 0.0022, -0.524/1e3},... % 273-363 K
                 'C2H4',  { 0.0037,  0.000    },... % 298     K
                 'C2H6',  { 0.0120, -0.601/1e3},... % 273-348 K
                 'C3H6',  { 0.0000,  0.000    },... % data not available
                 'C3H8',  { 0.0240, -0.702/1e3});   % 286-345 K

    % Model parameters to describe ionic dependence in the Sechenov relati-
    % on. The first parameter denotes the valence of the ion and the second
    % parameter is h_i [L/mol] and used to evaluate K [1/mol] (the Sechenov
    % constant) according to Eqn. (3) in the source:
    %
    %    K = SUM(n_i(h_i + h_G))
    %
    h_i = struct('H',    {+1, 0.0000},...
                 'Li',   {+1, 0.0754},...
                 'Na',   {+1, 0.1143},...
                 'K',    {+1, 0.0922},...
                 'Rb',   {+1, 0.0839},...
                 'Cs',   {+1, 0.0759},...
                 'NH4',  {+1, 0.0556},...
                 'Mg',   {+2, 0.1694},...
                 'OH',   {-1, 0.0839},...
                 'Cl',   {-1, 0.0318},...
                 'HCO3', {-1, 0.0967},...
                 'H2PO4',{-1, 0.0906},...
                 'HSO3', {-1, 0.0549},...
                 'CO3',  {-2, 0.1423},...
                 'HPO4', {-2, 0.1499},...
                 'SO3',  {-2, 0.1270},...
                 'SO4',  {-2, 0.1117},...
                 'PO4',  {-3, 0.2119});

    % CO2-specific temperature-corrected Sechenov relation
    h_G = h(1).(species) + h(2).(species).*(T - 298.15); % [L/mol]
    % Sechenov constant for the species in the electrolyte
    if strcmp('KHCO3',salt) == 1
        K = c.*((h_i(2).K + h_G) + (h_i(2).HCO3 + h_G))'; % [L/mol]
    elseif strcmp('KOH',salt) == 1
        K = c.*((h_i(2).K + h_G) + (h_i(2).OH + h_G))'; % [L/mol]
    elseif strcmp('KCl',salt) == 1
        K = c.*((h_i(2).K + h_G) + (h_i(2).Cl + h_G))'; % [L/mol]
    elseif strcmp('K2CO3',salt) == 1
        K = c.*(2*(h_i(2).K + h_G) + (h_i(2).CO3 + h_G))'; % [L/mol]
    % update Sechenov constant with simulated local concentrations
    % enter c as array with concentrations of HCO3(-), CO3(2-), OH(-), H(+) and K(+)
    elseif strcmp('update',salt) == 1
        K = c(1)*(h_i(2).HCO3 + h_G) + ...
            c(2)*(h_i(2).CO3 + h_G) + ...
            c(3)*(h_i(2).OH + h_G) + ...
            c(4)*(h_i(2).H + h_G) + ...
            c(5)*(h_i(2).K + h_G); % [L/mol]
    end
    % calculate c_sat_S in [mol/L]
    c_0     = repmat(henry(T,P,species)',1,length(c));
    c_sat_S = c_0./10.^(K);
end