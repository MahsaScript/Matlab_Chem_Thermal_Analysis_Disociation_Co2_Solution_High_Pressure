% Calculate the temperature- and ionic strength-dependent rate coefficient
% value of the CO2/carbonic acid/bicarbonate equilibrium. The system of di-
% ssociation reactions can be described by:
%
%                   k2f
%   H(+) + HCO3(-) <---> H2CO3
%          ʌ        k2b    ʌ
%      k3f | k3b       k1f | k1b
%          v               v
%              CO2 + H2O
%
% There is also the auxiliary reaction in alkaline medium according to:
%
%                   k4f
%          HCO3(-) <---> CO2 + OH(-)
%                   k4b
%
% The hydration rate of CO2 involves a combination of k1f and k3f under the
% assumption of instantaneous equilibrium of reaction #2. It is here retur-
% ned as combined parameter k1 = k1f + k3f. Furthermore, k2 = k4b*Kw as de-
% scribed in the appendix. The range of applicability is:
% 
%    278.15-308.15 K (determined at S = 33.77)
%    3.4-37.1 salinity units (determined at T = 298.15 K). 
%
% Salinity can be roughly equated with the mass of dissolved salts in grams
% per kilogram of seawater. Note that salinity values in the open oceans at
% midlatitudes typically fall between 34 and 36.

%  @T       temperature in [K]
%  @N       molarity in [M]
%  @type    string type of salt (i.e. 'KHCO3', 'KOH', or 'KCl')
%  @returns k1 - hydration rate coefficient of CO2 in [/s]
%           k2 - alkaline dissolution rate of CO2 in [L/mol/s]
function [k1,k2] = carbonatekin(T,c,type)
    % ionic strength and salinity of sea water
    I_ref = 0.7126; % [M]
    S_ref = 35;     % [g/kg]

    % calculate salinity in [g/kg]
    S = 0;
    if contains(type,'KHCO3') == 1 || contains(type,'KOH') == 1 || contains(type,'KCl') == 1
        S = S_ref/I_ref*c;
    elseif contains(type,'K2CO3') == 1
        S = S_ref/I_ref*3*c;
    end

    % in case vectors of salinity and temperature were specified
    %   S-dependence across a row of the matrix
    %   T-dependence across a column of the matrix
    [S,T] = meshgrid(S,T);
    k1 = exp( 1246.98 - 6.19e4./T - 183.0*log(T));
    k2 = exp( -930.13 + 3.10e4./T + 140.9*log(T) + 0.110*sqrt(S));
end