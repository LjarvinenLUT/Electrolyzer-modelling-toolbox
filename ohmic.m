% Ohmic overpotetial
% Inputs:   type - Electrolysis type, "PEM" or "alkaline"
%           conductivityModel - Used conductivity model 1-3
%           resistanceModel - Used resistance model
%           delta - Alkaline: Electrolyte layer thickness, PEM: Membrane thickness

function Uohm = ohmic(varargin)

    defaultConductivityModel = 1;
    defaultResistanceModel = 1;
    defaultType = 'Pem';

    parser = inputParser;
    addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
    addParameter(parser,'conductivityModel',defaultConductivityModel,@(x) isnumeric(x)&&isscalar(x));
    addParameter(parser,'resistanceModel',defaultResistanceModel,@(x) isnumeric(x)&&isscalar(x));
    addParameter(parser,'delta',@(x) isnumeric(x));
    addParameter(parser,'lambda',@(x) isnumeric(x));
    addParameter(parser,'m',@(x) isnumeric(x));
    addParameter(parser,'w',@(x) isnumeric(x));
    addParameter(parser,'Temperature',@(x) isnumeric(x));
    
    parse(parser, varargin{:});
    
    type = string(lower(parser.Results.type));
    conductivityModel = parser.Results.conductivityModel;
    resistanceModel = parser.Results.resistanceModel;
    delta = parser.Results.delta;
    lambda = parser.Results.lambda;
    m = parser.Results.m;
    w = parser.Results.w;
    T = parser.Results.Temperature;
    
    fprintf('\nOhmic overpotential calculation properties:\n')
    fprintf('Electrolyzer: %s\n', type)
    fprintf('Conductivity model: %d\n', conductivityModel)
    fprintf('Resistance model: %d\n', resistanceModel)
    
    %% Error checking
    
    switch resistanceModel
        case 1
            
        case 2
            if strcmpi(type,"PEM") && conductivityModel == 1
                    if ~isnumeric(lambda)
                        error('lambda has to be set for PEM conductivity model 1')
                    elseif  lambda <= 1
                        error('Lambda has to be greater than 1 when using conductivity model 1')
                    end
            end
            if strcmpi(type,"Alkaline")
                if conductivityModel == 1
                    if ~isnumeric(m)
                        error('m has to be set for alkaline conductivity model 1')
                    end
                elseif conductivityModel == 2
                    if ~isnumeric(w)
                        error('w has to be set for alkaline conductivity model 2')
                    end
                end
            end
            if  ~isnumeric(delta)
                error('Delta has to be set when using resistance model 2')
            end
            if ~isnumeric(T)
                error('Temperature has to be set when using resistance model 2')
            end
    end

    
    %%
    
    % Conductivity equations for PEM and Alkaline
    %   lambda -  water content
    %   m - molar concentration
    %   w - mass concentration wt%
    % Output
    %   sigma - specific conductivity (S/cm)
    if strcmpi(type,"Alkaline") && resistanceModel ~= 1
        switch conductivityModel
            case 1 % (Alkaline) Gilliam et al. "A review of specific conductivities of potassium hydroxide solutions for various concentrations and temperatures", 2007
                A = -2.041;  B = -0.0028;  C = 0.005332;
                D = 207.2;  E = 0.001043; F = -0.0000003;
                sigma = A*m + B*m^2 + C*m*T + D*(m/T) + E*m^3 + F*m^2*T^2;
            case 2 % (Alkaline, KOH) See et al. "Temperature and Concentration Dependence of the Specific Conductivity of Concentrated Solutions of Potassium Hydroxide"
                K1 = 0.279844803; K2 = -0.009241294; K3 = -0.000149660371;
                K4 = -0.000905209551; K5 = 0.000114933252; K6 = 0.1765;
                K7 = 0.0696648518; K8 = -28.9815658;
                sigma = K1*(100*w) + K2*T + K3 * T^2 + K4 * (T*100*w) ...
                    + K5*(T^2*(100*w)^K6) + K7*(T/(100*w)) + K8*((100*w)/T);
        end
    elseif strcmpi(type,"PEM") && resistanceModel ~= 1
        switch conductivityModel
            case 1 % (PEM) https://doi.org/10.1149/1.2085971 Springer et al. "Polymer Electrolyte Fuel Cell Model"
                sigma = 0.005139 * lambda - 0.00326 * exp(1268 * (1/303 - 1/T));
        end
    end
    
    switch resistanceModel
        case 1 % https://doi.org/10.1016/j.ijhydene.2015.03.164 "Electrochemical performance modeling of a proton exchange membrane electrolyzer cell for hydrogen energy"
            Uohm = @(r, j) r.*j;
        case 2
            Ri = delta/sigma;
            Uohm = @(r, j) (r + Ri).*j;  
    end
end