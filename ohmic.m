% Ohmic overpotetial
% Inputs:   D - Parameters
%           model - Used model reference, numeric, from which article
%           conductivityModel - Used conductivity model

function Uohm = ohmic(D, model, conductivityModel)
    
    r = D(6);
    j = D(5);
    
    Uohm = r.*j;

    % Conductivity equations for PEM and Alkaline
    %   lambda -  water content
    %   m - molar concentration
    %   w - mass concentration wt%
    % Output
    %   sigma - specific conductivity (S/cm)
    switch conductivityModel
        case 1 % (PEM) https://doi.org/10.1149/1.2085971 Springer et al. "Polymer Electrolyte Fuel Cell Model"
            sigma = (0.005139 * lambda - 0.0326 * exp(1268 * (1/303 - 1/T)));
        case 2 % (Alkaline) Gilliam et al. "A review of specific conductivities of potassium hydroxide solutions for various concentrations and temperatures", 2007
            A = -2.041;  B = -0.0028;  C = 0.005332;
            D = 207.2;  E = 0.001043; F = -0.0000003;
            sigma = A*m + B*m^2 + C*m*T + D*(m/T) + E*m^3 + F*m^2*T^2;
        case 3 % (Alkaline) See et al. "Temperature and Concentration Dependence of the Specific Conductivity of Concentrated Solutions of Potassium Hydroxide"
            % Todo: Figure out how this was formed
            sigma = -2.96396 - 0.02371*T + 0.12269*w + (5.7*e - 5)*T^2 ...
            + 0.00173*w^2 + (4.7*e - 4)*w - (3.6*e - 8)*T^3 + (2.7*e - 6)*w^3 ...
            - (8.9*e - 6)*T*w^2 + (2.4*e - 7)*T^2*w;
    end
end