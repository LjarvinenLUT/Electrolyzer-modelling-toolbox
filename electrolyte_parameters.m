% Function for determining liquid electrolyte properties
% Outputs:  psvEl   - Saturated vapor pressure of the electrolyte
%           aH2OEl  - Water activity
% Inputs:   T       - Measured temperature
%           m       - Electrolyte molality, mol/kg of solvent
%           electrolyte - Name of the electrolyte
%          

function [psvEl,aH2OEl] = electrolyte_parameters(T,m,electrolyte,varargin)

defaultModel = 1; % Default experimental parameters

parser = inputParser;
addRequired(parser,'T',@(x) isnumeric(x));
addRequired(parser,'m',@(x) isnumeric(x));
addRequired(parser,'electrolyte',@(x) ischar(x)||isstring(x))
addOptional(parser,'model',@(x) isnumeric(x)&&isscalar(x))

parse(parser,T,m,electrolyte,varargin{:});

model = parser.Results.model;


%% Case structure of models

switch model
    case 1 % Experimental fit, Balej, J., "Water Vapour Partial Pressures and Water Activities in Potassium and Sodium Hydroxide Solutions Over Wide Concentration and Temperature Ranges", J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
        if electrolyte == "KOH"
            if any(m < 0|m > 18)
                error('Molality for KOH has to be between 0 and 18 mol kg^-1.')
            end
%             if any(T<273.15|T>474)
%                 error('Experimental electrolyte data defined for temperatures between 273 K and 474 K.')
%             end
            
            % Parameters for partial pressure
            a = -0.01508*m - 0.0016788*m.^2 + 2.25887e-5*m.^3;
            b = 1 - 0.0012062*m + 5.6024e-4*m.^2 - 7.8228e-6*m.^3;
            
            % Activity
            aH2OEl = 10.^(-0.02255.*m + 0.001434.*m.^2 + (1.38*m - 0.9254*m.^2)./T);
            
        elseif electrolyte == "NaOH"
            if any(m < 0|m > 25)
                error('Molality for NaOH has to be between 0 and 25 mol kg^-1.')
            end
%             if any(T<273|T>474)
%                 error('Experimental electrolyte data defined for temperatures between 273 K and 474 K.')
%             end
            
            % Parameters for partial pressure
            a = -0.010986*m - 1.461e-3*m.^2 + 2.03528e-5*m.^3;
            b = 1 - 1.34141e-3*m + 7.07241e-4*m.^2 - 9.5362e-6*m.^3;
            
            % Activity
            aH2OEl = 10.^(-0.01332*m + 0.002542.*m.^2 - 3.06e-5.*m.^3 + (1.5827.*m - 1.5669.*m.^2 + 0.021296.*m.^3)./T);
            % Fixed the equation from the publication
            
        else
            error('Only KOH and NaOH defined as alkaline electrolytes for the experimental data.')
        end
        
        % Pure water vapor pressure
        psv = water_vapor_pressure(T,'model',3);
        
        % Electrolyte vapor pressure
        psvEl = psv.^b.*10.^(a);

    case 2 % Fitted model, Ursúa, A. and Sanchis, P., "Static–dynamic modelling of the electrical behaviour ofa commercial advanced alkaline water electrolyser", Int. J. Hydrogen Energy, 37(24):18598–18614, 2012
        
        
        
end