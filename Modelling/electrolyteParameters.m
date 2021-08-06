function [psvEl,aH2OEl] = electrolyteParameters(T,m,electrolyte,varargin)
% ELECTROLYTEPARAMETERS Calculate saturated water vapor pressure and water
% activity for a water solution of KOH or NaOH.
%
%   [psvEl,aH2OEl] = ELECTROLYTEPARAMETERS(T,m,electrolyte) calculates values
%       using equations from Balej, J., "Water Vapour Partial Pressures and 
%       Water Activities in Potassium and Sodium Hydroxide Solutions Over 
%       Wide Concentration and Temperature Ranges", 
%       J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
%       - Inputs:   T -- Temperature in Kelvin
%                   m -- Solution molality, mol of solute/kg of solvent
%                   electrolyte -- Electrolyte type: 1 = KOH
%                                                    2 = NaOH
%       -Outputs:   psvEl -- Saturated water vapor pressure for the
%       solution in absolute bar
%                   aH2OEl -- Water activity for the solution
%
%   [psvEl,aH2OEl] = ELECTROLYTEPARAMETERS(_,'model',model)
%       allows changing the model for one of the two contained:
%           #1 -- Equations presented by Balej J. (default)
%           #2 -- Ideal solution approximation
%
% See also WATERVAPORPRESSURE

defaultModel = 1; % Default experimental parameters

parser = inputParser;
addRequired(parser,'T',@(x) isnumeric(x)); % K, Temperature
addRequired(parser,'m',@(x) isnumeric(x)); % mol/kg of solvent, Molality
addRequired(parser,'electrolyte',@(x) isnumeric(x)&&isscalar(x)); % Electrolyte, 1 = KOH, 2 = NaOH
addOptional(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))

parse(parser,T,m,electrolyte,varargin{:});

model = parser.Results.model;


%% Case structure of models

switch model
    case 1 % Fitted model, Balej, J., "Water Vapour Partial Pressures and Water Activities in Potassium and Sodium Hydroxide Solutions Over Wide Concentration and Temperature Ranges", J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
        if electrolyte == 1 % KOH
            if any(m < 0|m > 18)
                error('Molality for KOH has to be between 0 and 18 mol kg^-1.')
            end
            if any(T<273.15|T>474)
                error('Electrolyte data defined for temperatures between 273 K and 474 K.')
            end
            
            % Parameters for partial pressure
            a = -0.01508*m - 1.6788e-3*m.^2 + 2.25887e-5*m.^3;
            b = 1 - 1.2062e-3*m + 5.6024e-4*m.^2 - 7.8228e-6*m.^3;
            
            % Activity
            aH2OEl = 10.^(-0.02255.*m + 0.001434.*m.^2 + (1.38*m - 0.9254*m.^2)./T);
            
        elseif electrolyte == 2 % NaOH
            if any(m < 0|m > 25)
                error('Molality for NaOH has to be between 0 and 25 mol kg^-1.')
            end
            if any(T<273|T>474)
                error('Electrolyte data defined for temperatures between 273 K and 474 K.')
            end
            
            % Parameters for partial pressure
            a = -0.010986*m - 1.461e-3*m.^2 + 2.03528e-5*m.^3;
            b = 1 - 1.34141e-3*m + 7.07241e-4*m.^2 - 9.5362e-6*m.^3;
            
            % Activity
            aH2OEl = 10.^(-0.01332*m + 0.002542.*m.^2 - 3.06e-5.*m.^3 + (1.5827.*m - 1.5669.*m.^2 + 0.021296.*m.^3)./T);
            % Fixed the equation from the publication
            
        else
            error('Only KOH and NaOH defined as alkaline electrolytes.')
        end
        
        % Pure water vapor pressure [bara]
        psv = waterVaporPressure(T,'model',2);
        
        % Electrolyte vapor pressure [bara]
        psvEl = psv.^b.*10.^(a);

        
    case 2 % Ideal solution approximation
        
        MH2O = (2*1.008 + 16)*1e-3; % kg/mol, molar mass of water
        
        molfracSolute = m*MH2O; % Molar fraction of solute
        molfracH2O = 1-molfracSolute; % Molar fraction of H2O

        % Pure water vapor pressure [bara]
        psv = waterVaporPressure(T,'model',2);
        
        % Electrolyte vapor pressure [bara]
        psvEl = psv*molfracH2O;
        
        % Water activity
        aH2OEl = psvEl./psv;
        
end