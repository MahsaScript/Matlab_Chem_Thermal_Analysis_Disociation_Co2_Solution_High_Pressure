
% Self-ionization constant of water, pKw = -log(Kw), across 0-100 degC eva-
% luated at saturated water vapour pressure. This does not contain the dep-
% endency on pressure or ionic strength of electrolyte, of which it is less
% strong of a function than temperature (i.e. |ΔpKw| > 0.3 for ΔT = 10 °C
% versus |ΔpKw| < 0.001 for ΔP = 1 bar). The data is fitted by a polynomial 
% equation of order 2 according to tabulated data in the CRC on page 5-71.
% The constant is defined as:
%    
%    H2O <---> H(+) + OH(-)  with Kw = [H(+)][OH(-)]
%    
% Note that pKa + pKb = pKw.
%
% On the other hand, the formula reported by Millero includes also ionic
% strength dependence through salinity, though possibly the range of appli-
% cability is, like mentioned in the file carbonateeq.m:
%    
%    273.15-323.15 K
%    1-50 salinity units
%

%  @T       temperature in [K]
%  @c       molarity [M]
%  @type    string of salt identity (e.g. 'KHCO3')
%  @returns pKw in [-]   
function pKw = selfionization(T,c,type)
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
    Kw = exp(148.9802 - 13847.26./T - 23.6521.*log(T) + (-5.977 + 118.67./T + 1.0495.*log(T)).*sqrt(S) - 0.01615*S);
    pKw      = -log10(Kw);
end