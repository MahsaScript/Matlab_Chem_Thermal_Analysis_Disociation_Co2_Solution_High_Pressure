% Calculates the Bjerrum plot (the fraction of CO2/HCO3(-)/CO3(2-) of the
% Dissolved Organic Content (DIC) as a function of pH). This is calculated
% on the basis of specified equilibrium coefficients.

%  @Kw      water dissociation coefficient in [M2]
%  @Kb1     base equilibrium coefficient for the CO2/HCO3(-) pair in [M]
%  @Kb2     base equilibrium coefficient for the HCO3(-)/CO3(2-) pair in [M]
%  @returns vectors of pH in [-]
%           and the associated fraction of CO2 [-]
%                                          HCO3(-) [-]
%                                          CO3(2-) [-]
function [pH,f_CO2,f_iHCO3,f_iCO3] = bjerrum(Kw,Kb1,Kb2)
    pH      = linspace(0,14,150);
    pKw     = -log10(Kw);
    c_iOH   = 10.^(-pKw + pH);
    r_1     = Kb1./c_iOH;               % ratio of CO2 to HCO3(-)
    r_2     = Kb2./c_iOH;               % ratio of HCO3(-) to CO3(2-)
    r_3     = Kb1*Kb2./(c_iOH.^2);      % ratio of CO2 to CO3(2-)
    f_CO2   = 1./(1 + 1./r_1 + 1./r_3); % fraction CO2 to DIC
    f_iHCO3 = 1./(1 + r_1 + 1./r_2);    % fraction HCO3(-) to DIC
    f_iCO3  = 1./(1 + r_3 + r_2);       % fraction CO3(2-) to DIC

    % plot colours
    b = [0 0.4470 0.7410];
    o = [0.8500 0.3250 0.0980];
    y = [0.9290 0.6940 0.1250];

    fs = 13;
    semilogy(pH,f_CO2,'color',b,LineWidth = 2); hold on;
    semilogy(pH,f_iHCO3,'color',o,LineWidth = 2); hold on;
    semilogy(pH,f_iCO3,'color',y,LineWidth = 2); hold on;
    set(gca,'FontSize',fs);
    xlim([0 14]); ylim([1e-3 1.1]);
    xlabel('pH (-)','Fontsize',fs); ylabel('{\it c_i} /{\it c}_{DIC}','FontSize',fs);
    legend('{\it i} = CO_{2}','{\it i} = HCO_{3}^{-}','{\it i} = CO_{3}^{2-}','location','northwest');
end