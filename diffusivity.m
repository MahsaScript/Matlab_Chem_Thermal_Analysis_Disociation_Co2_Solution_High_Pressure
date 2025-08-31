% Diffusion coefficients of chemical species corrected for temperature and
% viscosity of the ionic solution based on the Stokes-Einstein equation:
%
%    D_0*u_0/T_0 = D*u/T
%
% The diffusivities of ions at infinite dilution can be found in the source
% on pages 5-77 to 5-79. The diffusivities of gases in water can be found
% on page 6-250. The diffusivities of gases in excess air can be found on
% page 6-248.
%
% Diffusivities are modified according to the porosity of a random porous
% medium made up of spherical particles. The Bruggeman correlation is invo-
% ked for this:
%
%    D_eff = D*e^(3/2)
%
% Ion mobilities are found according to the Einstein-Smoluchowski equation
%
%    u_i = D_i*z_i*F/(R*T)
%

%  @T       vector or scalar of temperature in [K]
%  @u       vector of scalar viscosity of electrolyte in [Pa s]
%  @u       vector of scalar viscosity of air in [Pa s]
%  @e       scalar of porosity (e = 0 if non-porous) in [-]
%  @returns structure of diffusion coefficients in water
%                        diffusion coefficients in air
%                        ion valence
%                        ionic mobility
function D_0 = diffusivity(T,u,u_air,e)

    % column 1: diffusion coefficients at infinite dilution in water [m2/s] at 298.15 K
    % column 2: diffusion coefficients at infinite dilution in air [m2/s] at 293.15 K
    % column 3: valence of ion
    % column 4: ion mobility in m2/s/V
    % field names preceded by 'i' indicate ionic species
    D_0 = struct('CO2',  {1.91e-9, 0.160e-4, 0,0},...
                 'CO',   {1e-9    ,0.208e-4, 0,0},...
                 'H2',   {5.11e-9 ,0.756e-4, 0,0},...
                 'NH3',  {1.5e-9  ,0.280e-4, 0,0},... % at 293.15
                 'iNH4', {1.957e-9,0       , 1,0},...
                 'iH',   {9.311e-9,0       , 1,0},...
                 'iOH',  {5.273e-9,0       ,-1,0},...
                 'iHCO3',{1.185e-9,0       ,-1,0},...
                 'iCO3' ,{0.923e-9,0       ,-2,0},...
                 'iK'   ,{1.957e-9,0       , 1,0},...
                 'iHCOO',{1.454e-9,0       ,-1,0},...
                 'iNO3', {1.902e-9,0       ,-1,0},...
                 'iCl',  {2.032e-9,0       ,-1,0});     

    % viscosity of pure H2O at 298.15 K and 1 bara
    u_water_0 = 8.90e-4; % [Pa s]
    % viscosity of air at 293.15 K and 1 bara
    u_air_0   = 1.82e-5; % [Pa s]
    % overwrite the previous structure with corrected diffusivities
    fields = fieldnames(D_0);
    for i  = 1:numel(fields)
        if strcmp(fields{i},'NH3') == 1
            D_0(1).(fields{i}) = D_0(1).(fields{i})*T/293.15*u_air_0./u.*e.^(3/2);
        else
            D_0(1).(fields{i}) = D_0(1).(fields{i})*T/298.15*u_water_0./u.*e.^(3/2);
            D_0(2).(fields{i}) = D_0(2).(fields{i})*T/293.15*u_air_0./u_air.*e.^(3/2);
        end
        D_0(4).(fields{i}) = D_0(1).(fields{i})*abs(D_0(3).(fields{i}))*96485./(8.3145.*T);
    end
end