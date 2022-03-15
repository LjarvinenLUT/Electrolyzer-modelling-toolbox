function Uohm = ohmic(varargin)
% OHMIC  Create a func object for calculation of ohmic overpotential
% for modelling of water electrolysis.
%
%   Uohm = OHMIC() creates a func object with simple Ohm's law, resistance
%           model #1.
% 
%   Uohm = OHMIC('resistanceModel',rm) uses resistance model defined by rm.
%           Available models are:
%               #1 -- Simple Ohm's law with total resistance. Output func
%                       object contains the following fields in its
%                       workspace:
%                       Variables: current -- Current density
%                       Parameters: r -- Area specific resistance
%               #2 -- Ohm's law with electrolyte specific resistance
%               separated.
%                   - Requires electrolysis type as a name-value pair:
%
%                     Uohm = OHMIC(_,'type',e)
%                     
%                     Available types are 'PEM' and 'alkaline'.
%
%                   - Conductivity model for alkali electrolysis can be
%                     changed with name-value pair:
%
%                     Uohm = OHMIC(_,'conductivityModel',cm)
%                       
%                     Available models are #1 and #2.
%
%                   - Parameters required for resistance model #2:
%                       delta -- electrolyte/membrane thickness
%                       lambda -- PEM membrane water content
%                       m -- Alkaline solution molality or molarity (TODO: check)
%                       w -- Alkaline solution mass concentration
%                       T -- Temperature
% 
%                   - Output func objects workspace contains the following 
%                       fields in workspace:
%                       Variables: 
%                           current -- Current density
%                           R_ionic -- Ionic resistance
%                       Parameters: 
%                           r_electronics -- Area specific resistance of
%                               the electronics only
%                      
% 
%   See also NERNST, ACTIVATION, CONCENTRATION, FUNC
    
    defaultConductivityModel = 1;
    defaultResistanceModel = 1;
    defaultType = 'pem';

    parser = inputParser;
    addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
    addParameter(parser,'conductivityModel',defaultConductivityModel,@(x) isnumeric(x)&&isscalar(x));
    addParameter(parser,'resistanceModel',defaultResistanceModel,@(x) isnumeric(x)&&isscalar(x));
    
    parse(parser, varargin{:});
    
    type = string(lower(parser.Results.type));
    conductivityModel = parser.Results.conductivityModel;
    resistanceModel = parser.Results.resistanceModel;
    Workspace.Variables = struct('current',[],'T',[]);

    
    fprintf('\nOhmic overpotential modelling properties:\n')
    
    %%
    
    

    
    switch resistanceModel
        case 1 % https://doi.org/10.1016/j.ijhydene.2015.03.164 "Electrochemical performance modeling of a proton exchange membrane electrolyzer cell for hydrogen energy"
            
            resistanceModelStr = resistanceModel + " -- Total cell resistance, combined electronic and ionic components";
            fprintf('Resistance model: %s\n', resistanceModelStr)
            
            Workspace.Parameters.r = [];
            Fitlims.r = {0,1,1000};
            funcHandle = @(Workspace) Workspace.Parameters.r.*Workspace.Variables.current;
            Workspace.Dependencies.warning_r = "if ~isempty(Workspace.Parameters.r)&&changed.T;"+...
                                "warning('Temperature dependency of the area specific resistivity (r) is not taken into account in the models. Changing the temperature without a new fit will result in unrealistic behaviour.');"+...
                                "end;";
        case 2
            % Conductivity equations for PEM and Alkaline
            %   lambda -  membrane water content (for PEM)
            %   m - molar concentration (for alkaline)
            %   delta -- electrolyte thickness (TODO: in which units)
            % Output
            %   sigma - specific conductivity (S/cm)
            
            % Printing explanation
            resistanceModelStr = resistanceModel + " -- Separated electronic and ionic resistance components";
            fprintf('Resistance model: %s\n', resistanceModelStr)
            fprintf('Electrolyzer type: %s\n', type)
            fprintf('TODO: Check units for ohmic overpotential parameters!\n')
            if strcmpi(type,"Alkaline")
                fprintf('Conductivity model: %d\n', conductivityModel)
            end
            
            Workspace.Variables.sigma = [];
            Workspace.Variables.delta = [];  
            
            switch type
                case "alkaline"
                    warning('Ohmic overpotential using resistance model 2 in alkaline is currently only available for KOH electrolyte!')
                    switch conductivityModel
                        case 1 % (Alkaline) Gilliam et al. "A review of specific conductivities of potassium hydroxide solutions for various concentrations and temperatures", 2007
%                             A = -2.041;  B = -0.0028;  C = 0.005332;
%                             D = 207.2;  E = 0.001043; F = -0.0000003;
%                             sigma = A*M + B*M^2 + C*M*T + D*(M/T) + E*M^3 + F*M^2*T^2;
                            Workspace.Dependencies.sigma = "Workspace.Variables.sigma = "...
                                + "-2.041*Workspace.Variables.Molarity + "...
                                + "-0.0028*Workspace.Variables.Molarity.^2 + "...
                                + "0.005332*Workspace.Variables.Molarity.*Workspace.Variables.T + "...
                                + "207.2*(Workspace.Variables.Molarity./Workspace.Variables.T) + "...
                                + "0.001043*Workspace.Variables.Molarity.^3 + "...
                                + "-0.0000003*Workspace.Variables.Molarity.^2*Workspace.Variables.T.^2;"; % Dependency format
                        case 2 % (Alkaline, KOH) See et al. "Temperature and Concentration Dependence of the Specific Conductivity of Concentrated Solutions of Potassium Hydroxide"
%                             K1 = 0.279844803; K2 = -0.009241294; K3 = -0.000149660371;
%                             K4 = -0.000905209551; K5 = 0.000114933252; K6 = 0.1765;
%                             K7 = 0.0696648518; K8 = -28.9815658;
%                             sigma = K1*(100*w) + K2*Workspace.Variables.T + K3 * Workspace.Variables.T^2 + K4 * (Workspace.Variables.T*100*w) ...
%                                 + K5*(Workspace.Variables.T^2*(100*w)^K6) + K7*(Workspace.Variables.T/(100*w)) + K8*((100*w)/Workspace.Variables.T);
                            Workspace.Dependencies.sigma =  "Workspace.Variables.sigma = "...
                                + "0.279844803*(100*Workspace.Variables.wtfrac) + "...
                                + "-0.009241294*Workspace.Variables.T + "...
                                + "-0.000149660371*Workspace.Variables.T.^2 + "...
                                + "-0.000905209551*(Workspace.Variables.T.*(100.*Workspace.Variables.wtfrac)) +"...
                                + "0.000114933252*(Workspace.Variables.T.^2*(100*Workspace.Variables.wtfrac).^0.1765) +"...
                                + "0.0696648518*(Workspace.Variables.T./(100*Workspace.Variables.wtfrac)) +"...
                                + "-28.9815658*((100*Workspace.Variables.wtfrac)./Workspace.Variables.T);"; % Dependency format
                    end
                case "pem"
                    % (PEM) https://doi.org/10.1149/1.2085971 Springer et al. "Polymer Electrolyte Fuel Cell Model"
%                     sigma = (0.005139 * lambda - 0.00326) * exp(1268 * (1/303 - 1./Workspace.Variables.T));
                    Workspace.Variables.lambda = [];
                    Workspace.Dependencies.sigma =  "Workspace.Variables.sigma = "...
                        + "(0.005139*Workspace.Variables.lambda - 0.00326)"...
                        + "*exp(1268*(1/303-1./Workspace.Variables.T));"; % Dependecy format
            end
            
                      
            Workspace.Dependencies.R_ionic = "Workspace.Variables.R_ionic = Workspace.Variables.delta./Workspace.Variables.sigma;";
            Workspace.Parameters.r_electronics = [];
            Fitlims.r_electronics = {0,1,1000};
            funcHandle = @(Workspace) (Workspace.Parameters.r_electronics + Workspace.Variables.R_ionic).*Workspace.Variables.current;
    end
    
    Workspace.Constants = struct();
    
    Uohm = func(funcHandle,Workspace,Fitlims);
    
end