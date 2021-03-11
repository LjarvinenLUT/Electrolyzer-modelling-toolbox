% Ohmic overpotetial
% Inputs:   delta - Alkaline: Electrolyte layer thickness, PEM: Membrane thickness
%           type - Electrolysis type, "PEM" or "alkaline"
%           conductivityModel - Used conductivity model 1-3
%           resistanceModel - Used resistance model

function Uohm = ohmic(delta, type, conductivityModel, resistanceModel) 

    parser = inputParser;
    addRequired(parser,'r',@(x) isnumeric(x));
    addRequired(parser,'j',@(x) isnumeric(x));
    addRequired(parser,'d',@(x) isnumeric(x));
    addRequired(parser,'type',@(x) ischar(x)||isstring(x));
    addRequired(parser,'conductivityModel',@(x) isnumeric(x)&&isscalar(x));
    addRequired(parser,'resistanceModel',@(x) ischar(x)||isstring(x));
    
    parse(parser,r, j, d, type, conductivityModel, resistanceModel);
    
    %% Error checking
    if strcmpi(type,"Alkaline") && conductivityModel == 1
        error('Conductivity model 1 is only for PEM electrolysers')
    elseif strcmpi(type,"PEM") && conductivityModel ~= 1
        error('Use conduction model 1 for PEM electrolyser')    
    end
    
    %%
    
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
        case 3 % (Alkaline, KOH) See et al. "Temperature and Concentration Dependence of the Specific Conductivity of Concentrated Solutions of Potassium Hydroxide"
            K1 = 0.279844803; K2 = -0.009241294; K3 = -0.000149660371;
            K4 = -0.000905209551; K5 = 0.000114933252; K6 = 0.1765;
            K7 = 0.0696648518; K8 = -28.9815658;
            sigma = K1*(100*w) + K2*T + K3 * T^2 + K4 * (T*100*w) ...
                + K5*(T^2*(100*w)^K6) + K7*(T/(100*w)) + K8*((100*w)/T);
%         case 3 % (Alkaline) See et al. "Temperature and Concentration Dependence of the Specific Conductivity of Concentrated Solutions of Potassium Hydroxide"
%             % Todo: Figure out how this was formed
%             sigma = -2.96396 - 0.02371*T + 0.12269*w + (5.7*e - 5)*T^2 ...
%             + 0.00173*w^2 + (4.7*e - 4)*w - (3.6*e - 8)*T^3 + (2.7*e - 6)*w^3 ...
%             - (8.9*e - 6)*T*w^2 + (2.4*e - 7)*T^2*w;
    end
    
    switch resistanceModel
        case 1 % https://doi.org/10.1016/j.ijhydene.2015.03.164 "Electrochemical performance modeling of a proton exchange membrane electrolyzer cell for hydrogen energy"
            Uohm = @(r, j) r.*j;
        case 2
            Ri = delta/sigma;
            Uohm = @(r, j) (r + Ri).*j;
        case 3     
            %Uohm = @(r, j) (Rp_a + Re_a + Rin_a).*j*A + (Rp_c + Re_c + Rin_c).*j*A + d_m.*(j)./(sigma);
    end
end