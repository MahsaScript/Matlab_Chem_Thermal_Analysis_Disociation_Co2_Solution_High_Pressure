% Calculate the temperature- and ionic strength-dependent pKa values of the
% carbonate equilibria. The CO2/HCO3(-) couple is described by index 1, and
% the HCO3(-)/CO3(2-) couple by index 2 as below:
% 
%    Ka1 = [H(+)][HCO3(-)]/[CO2]     and pKa1 = -log10(Ka1)
%    Ka2 = [H(+)][CO3(2-)]/[HCO3(-)] and pKa2 = -log10(Ka2)
% 
% Reported range of applicability is:
%    
%    273.15-323.15 K
%    1-50 salinity units
%

%  @T       vector of temperature in [K]
%  @c       vector of salt molarity in [M]
%  @type    string type of salt (i.e. 'KHCO3' or 'KOH')
%  @returns Ka values in [mol/L]
%           pKa1_0 and pKa2_0 are temperature-corrected (vector)
%           pKa1   and pKa2   are temperature- and salinity-corrected (matrix)
function [Ka1,Ka2,Ka1_0,Ka2_0] = carbonateeq(T,c,type)
    % temperature-dependent values of pKa in pure water
    pKa1_0 = (-126.34048 + 6320.813./T + 19.568224*log(T))'; Ka1_0 = 10.^(-pKa1_0);
    pKa2_0 = (-90.18333 + 5143.692./T + 14.613358*log(T))';  Ka2_0 = 10.^(-pKa2_0);

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
    [S,T] = meshgrid(S,T);
    % adjustable parameters 
    A1 = 13.4191*sqrt(S) + 0.0331*S - 5.33e-5.*S.^2;
    A2 = 21.0894*sqrt(S) + 0.1248*S - 3.687e-4.*S.^2;
    B1 = -530.123*sqrt(S) - 6.103*S;
    B2 = -772.483*sqrt(S) - 20.051*S;
    C1 = -2.06950*sqrt(S);
    C2 = -3.3336*sqrt(S);
    % in case vectors of salinity and temperature were specified
    %   S-dependence across a row of the matrix
    %   T-dependence across a column of the matrix
    pKa1_0 = repmat(pKa1_0,1,length(c));
    pKa2_0 = repmat(pKa2_0,1,length(c));
    pKa1 = pKa1_0 + A1 + B1./T + C1.*log(T); Ka1 = 10.^(-pKa1); 
    pKa2 = pKa2_0 + A2 + B2./T + C2.*log(T); Ka2 = 10.^(-pKa2);
end