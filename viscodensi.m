
% Returns the  viscosity and density of salt solutions at 293.15 K as func-
% tion of molarity. The data from the source is fitted by a polynomial equ-
% ation of order 2 according to tabulated data in the source, pages 5-123 
% to 5-148. The script provides support for the following salt solutions:
%  - potassium bicarbonate (string = 'KHCO3')
%  - potassium carbonate (string = 'K2CO3')
%  - potassium hydroxide (string = 'KOH')
%  - potassium chloride (string = 'KCl')
%  - potassium sulfate (string = 'K2SO4')

%  @T       scalar or vector of temperature in [K]
%  @c       scalar or vector of molarity in [M]
%  @salt    string element
%  @returns scalar or vector of viscosity in [Pa s]
%                               density in [kg/m3]
function [u,p,u_air,p_air] = viscodensi(T,c,salt)
    % viscosity and density of water at sat. vapour pressure and 293.15 K
    T_0    = 293.15; % [K]
    u_0    = 1.002;  % [mPa s]
    p_0    = 0.998;  % [g/mL]
    % molarity in [mol/L]
    c_data = struct('KHCO3', {[0.000 0.050 0.100 0.202 0.305 0.409 0.515 0.622 0.730 0.840 0.951 1.064 1.293 1.528 1.770 2.017 2.272 2.533 2.801]},...
                    'K2CO3', {[0.000 0.036 0.073 0.147 0.223 0.299 0.378 0.457 0.538 0.620 0.704 0.789 0.963 1.144 1.330 1.523 1.722 2.139 2.584 3.057 3.559 4.093 5.573]},...
                    'KOH',   {[0.000 0.089 0.179 0.362 0.548 0.736 0.929 1.124 1.322 1.524 1.729 1.938 2.365 2.806 3.261 3.730 4.213 4.711 5.223 5.750 6.293 6.851 9.896 13.389]},...
                    'KCl',   {[0.000 0.067 0.135 0.271 0.409 0.549 0.691 0.835 0.980 1.127 1.276 1.426 1.733 2.048 2.370 2.701 3.039 3.386 3.742]},...
                    'K2SO4', {[0.000 0.037 0.058 0.116 0.176 0.237 0.298 0.360 0.424 0.488 0.554 0.620]});
    % viscosity in [Pa s]
    u_data = struct('KHCO3', {[u_0   1.009 1.015 1.027 1.040 1.053 1.067 1.081 1.096 1.112 1.128 1.145 1.183 1.224 1.270 1.319 1.373 1.432 1.497]/1e3},...
                    'K2CO3', {[u_0   1.013 1.025 1.048 1.071 1.094 1.119 1.146 1.174 1.204 1.235 1.269 1.339 1.414 1.497 1.594 1.707 1.978 2.331 2.834 3.503 4.360  9.369]/1e3},...
                    'KOH',   {[u_0   1.010 1.019 1.038 1.058 1.079 1.102 1.126 1.151 1.177 1.205 1.233 1.294 1.361 1.436 1.521 1.619 1.732 1.861 2.006 2.170 2.357  3.879  7.892]/1e3},...
                    'KCl',   {[u_0   1.000 0.999 0.999 0.998 0.997 0.996 0.994 0.992 0.990 0.989 0.988 0.990 0.994 0.999 1.004 1.012 1.024 1.040]/1e3},...
                    'K2SO4', {[u_0   1.006 1.011 1.021 1.033 1.045 1.058 1.072 1.087 1.102 1.117 1.132]/1e3});
    p_data = struct('KHCO3', {[p_0   1.001 1.005 1.011 1.018 1.025 1.031 1.038 1.045 1.051 1.058 1.065 1.079 1.093 1.107 1.122 1.137 1.157 1.169]*1e3},...
                    'K2CO3', {[p_0   1.003 1.007 1.016 1.025 1.034 1.044 1.053 1.062 1.072 1.081 1.090 1.109 1.129 1.149 1.169 1.190 1.232 1.276 1.320 1.367 1.412 1.540]*1e3},...
                    'KOH',   {[p_0   1.003 1.007 1.016 1.024 1.033 1.042 1.051 1.060 1.069 1.078 1.087 1.106 1.125 1.144 1.163 1.182 1.201 1.221 1.241 1.261 1.281 1.388 1.502]*1e3},...
                    'KCl',   {[p_0   1.001 1.005 1.011 1.017 1.024 1.030 1.037 1.043 1.050 1.057 1.063 1.077 1.091 1.104 1.118 1.133 1.147 1.162]*1e3},...
                    'K2SO4', {[p_0   1.002 1.006 1.014 1.022 1.031 1.039 1.047 1.055 1.064 1.072 1.081]*1e3});
    
    % temperature-dependence of viscosity and density of water in [K], [Pa s], and [kg/m3] respectively
    T_data      = [273.15 283.15 T_0 298.15 303.15 313.15 323.15 333.15 343.15 353.15 363.15 373.15];
    uT_data     = [1.791 1.306 u_0 0.890 0.797 0.653 0.547 0.466 0.404 0.354 0.314 0.283]/1e3;
    pT_data     = [1.000 1.000 p_0 0.997 0.996 0.992 0.988 0.983 0.978 0.972 0.965 0.959]*1e3;
    % temperature-dependence of viscosity and density of air in [K], [Pa s], and [kg/m3] respectively
    T_data_air  = [293.15 303.15 313.15 333.15 353.15 373.15];
    uT_data_air = [  1.82   1.87   1.92   2.01   2.10   2.18]/1e5;
    pT_data_air = [  1.205  1.165  1.127  1.060  1.000  0.946];

    % fitting functions for viscosity and density as function of temperature for the electrolyte
    poly_u  = polyfitB(c_data(1).(salt), u_data(1).(salt), 3, u_0*1e-3);
    u       = poly_u(1)*c.^3 + poly_u(2)*c.^2 + poly_u(3)*c + poly_u(4);
    poly_p  = polyfitB(c_data(1).(salt), p_data(1).(salt), 3, p_0*1e3);
    p       = poly_p(1)*c.^3 + poly_p(2)*c.^2 + poly_p(3)*c + poly_p(4);

    % additional correction for temperature (assumes viscosity and density scale similarly as for water)
    poly_uT = polyfitB(T_data, uT_data, 5, uT_data(1));
    u       = u*(poly_uT(1)*T.^5 + poly_uT(2)*T.^4 + poly_uT(3)*T.^3 + poly_uT(4)*T.^2 + poly_uT(5)*T + poly_uT(6))/...
                (poly_uT(1)*T_0^5 + poly_uT(2)*T_0^4 + poly_uT(3)*T_0^3 + poly_uT(4)*T_0.^2 + poly_uT(5)*T_0 + poly_uT(6));
    poly_pT = polyfitB(T_data, pT_data, 3, uT_data(1));
    p       = p*(poly_pT(1)*T.^3 + poly_pT(2)*T.^2 + poly_pT(3)*T + poly_pT(4))/...
                (poly_pT(1)*T_0^3 + poly_pT(2)*T_0^2 + poly_pT(3)*T_0 + poly_pT(4));

    % fitting functions for viscosity and density as function of temperature for moisture-free air
    poly_uT_air = polyfit(T_data_air, uT_data_air, 1);
    u_air       = poly_uT_air(1)*T + poly_uT_air(2);
    poly_pT_air = polyfit(T_data_air, pT_data_air, 1);
    p_air       = poly_pT_air(1)*T + poly_pT_air(2);

%     figure;
%     plot(T_data, pT_data,'*k'); hold on;
%     plot(T_data,polyval(poly_pT,T_data));
end