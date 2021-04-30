% Pressure correction for Nernst equation
% Inputs:   T - Measured temperature
%           p1 - Parameter 1: 
%               for "PEM" cathode pressure
%               for "alkaline" system pressure
%           p2 - Parameter 2: 
%               for "PEM" anode pressure
%               for "alkaline" electrolyte molality
%           type - Electrolysis type, "PEM" or "alkaline"
%           electrolyte - Electrolyte used for "alkaline"

function Ucor = nernst_pressure_correction(vars,varargin)
    
    defaultElectrolyte = 'KOH';

    parser = inputParser;
    addRequired(parser,'vars',@(x) isstruct(x));
    addParameter(parser,'type',@(x) ischar(x)||isstring(x));
    addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
    
    parse(parser,vars,varargin{:});
    
    type = string(lower(parser.Results.type));
    electrolyte = string(parser.Results.electrolyte);
    
    %% Errors
    if strcmp(type,"alkaline")
        if ~strcmp(electrolyte,"KOH")&&~strcmp(electrolyte,"NaOH")
            error('Only KOH and NaOH defined as alkaline electrolytes')
        end
    elseif ~strcmp(type,"pem")
        error('Only PEM and alkaline electrolysis defined for Nernst equation.')
    end
    %%
    [F,R,n_e] = get_constants;

    %%

    switch type
        case "pem"
            psv = water_vapor_pressure(vars.T); % Vapor pressure from Antoine equation
            
            pH2 = vars.p1 - psv; % bar, Hydrogen partial pressure
            pO2 = vars.p2 - psv; % bar, Oxygen partial pressure
            
            if any(pH2 < 0 | pO2 < 0)
                error('Pressure too low! Anode or cathode pressure lower than saturated vapor pressure.')
            end
            
            % Nernst equation pressure correction
            Ucor = (R.*vars.T)/(n_e*F).*log(pH2.*pO2.^(1/2));
            
        case "alkaline"
            
            [psvEl,aH2OEl] = electrolyte_parameters(vars.T,vars.m,electrolyte);
            
            ps = vars.p - psvEl; % bar, Partial pressure of hydrogen and oxygen
            
            if any(ps < 0)
                error('Pressure too low! System pressure lower than saturated vapor pressure of the electrolyte solution.')
            end
            
            % Nernst equation pressure correction
            Ucor = (R.*vars.T)/(n_e*F).*log(ps.^(3/2)./aH2OEl);
    end
end