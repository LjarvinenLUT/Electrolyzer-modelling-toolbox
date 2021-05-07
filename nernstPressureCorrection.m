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

function Ucor = nernstPressureCorrection(T,p1,p2,varargin)
    
    defaultElectrolyte = 'KOH';

    parser = inputParser;
    addRequired(parser,'T',@(x) isnumeric(x));
    addRequired(parser,'p1',@(x) isnumeric(x));
    addRequired(parser,'p2',@(x) isnumeric(x));
    addParameter(parser,'type',@(x) ischar(x)||isstring(x));
    addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
    
    parse(parser,T,p1,p2,varargin{:});
    
    Variables.T = T;
    type = string(lower(parser.Results.type));
    if strcmp(type,"alkali")
        type = "alkaline";
    end
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
    [Constants.F,Constants.R,Constants.n_e] = get_constants;

    %%

    switch type
        case "pem"
            psv = water_vapor_pressure(T); % Vapor pressure from Antoine equation
            
            Variables.pH2 = p1 - psv; % bar, Hydrogen partial pressure
            Variables.pO2 = p2 - psv; % bar, Oxygen partial pressure
            
            if any(Variables.pH2 < 0 | Variables.pO2 < 0)
                error('Pressure too low! Anode or cathode pressure lower than saturated vapor pressure.')
            end
            
            Workspace = struct('Constants',Constants,'Variables',Variables);
            
            % Nernst equation pressure correction
            funcHandle = @(Workspace) (Workspace.Constants.R.*Workspace.Variables.T)/(Workspace.Constants.n_e*Workspace.Constants.F).*log(Workspace.Variables.pH2.*Workspace.Variables.pO2.^(1/2));
            
        case "alkaline"
            
            [psvEl,Variables.aH2OEl] = electrolyte_parameters(Variables.T,p2,electrolyte);
            
            Variables.ps = p1 - psvEl; % bar, Partial pressure of both hydrogen and oxygen in the system
            
            if any(Variables.ps < 0)
                error('Pressure too low! System pressure lower than saturated vapor pressure of the electrolyte solution.')
            end
            
            Workspace = struct('Constants',Constants,'Variables',Variables);
            
            % Nernst equation pressure correction
            funcHandle = @(Workspace) (Workspace.Constants.R.*Workspace.Variables.T)/(Workspace.Constants.n_e*Workspace.Constants.F).*log(Workspace.Variables.ps.^(3/2)./Workspace.Variables.aH2OEl);
    end
    
    Ucor = func(funcHandle,Workspace);
end