function Ucor = nernstPressureCorrection(T,p1,p2,varargin)
% NERNSTPRESSURECORRECTION  Create a func object for pressure correction 
% term of Nernst equation in modelling of water electrolysis.
%
%   Ucor = NERNSTPRESSURECORRECTION(T,p1,p2,'type',t) creates func object
%           for electrolyzer type defined by variable t. 
%           Electrolyzer types availabel are:
%               PEM -- Polymer electrolyte membrane electrolysis
%               Alkaline -- Alkaline electrolysis.
%           Inputs: T -- Measured temperature
%                   p1 -- Parameter 1: 
%                       for PEM: cathode pressure, in bar
%                       for alkaline: system pressure, in bar
%                   p2 - Parameter 2: 
%                       for PEM: anode pressure, in bar
%                       for alkaline: electrolyte molality, in mol/kg of solvent
% 
%   Ucor = NERNSTPRESSURECORRECTION(_,'type','alkaline','electrolyte',e) changes the 
%           electrolyte when alkaline electrolysis is considered.
%           Available electrolytes are:
%               KOH -- Potassium hydroxide (default)
%               NaOH -- Sodium hydroxide
%   
%   The output func object has the following fields in workspace
%   structures, depending on the type of electrolyzer:
%       Constants (obtained by the function automatically):
%           F -- Faraday's constant
%           R -- Universal gas constant
%           n_e -- Number of electrons transferred in a single
%                   electrochemical water splitting reaction
%       Varaibles:
%           T -- Temperature in kelvins
%           pO2 -- Oxygen partial pressure (for PEM only)
%           pH2 -- Hydrogen partial pressure (for PEM only)
%           ps -- System pressure (for alkaline only)
%           aH2OEl -- Water activity in electrolyte solution
% 
%   See also NERNST, REVERSIBLE, FUNC, WATERVAPORPRESSURE,
%   ELECTROLYTEPARAMETERS
    
%% Parse input
    defaultElectrolyte = 'KOH';

    parser = inputParser;
    addRequired(parser,'T',@(x) isnumeric(x));
    addRequired(parser,'p1',@(x) isnumeric(x));
    addRequired(parser,'p2',@(x) isnumeric(x));
    addParameter(parser,'type',@(x) ischar(x)||isstring(x));
    addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
    
    parse(parser,T,p1,p2,varargin{:});
    
    Workspace.Variables = struct('T',T);
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
    %% Constants
    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = getConstants;

    %% Define the func object

    switch type
        case "pem"
            psv = waterVaporPressure(Workspace.Variables.T); % Vapor pressure from Antoine equation
            
            Workspace.Variables.pH2 = p1 - psv; % bar, Hydrogen partial pressure
            Workspace.Variables.pO2 = p2 - psv; % bar, Oxygen partial pressure
            
            if any(Workspace.Variables.pH2 < 0 | Workspace.Variables.pO2 < 0)
                error('Pressure too low! Anode or cathode pressure lower than saturated vapor pressure.')
            end
                      
            % Nernst equation pressure correction
            funcHandle = @(Workspace) (Workspace.Constants.R.*Workspace.Variables.T)/(Workspace.Constants.n_e*Workspace.Constants.F).*log(Workspace.Variables.pH2.*Workspace.Variables.pO2.^(1/2));
            
        case "alkaline"
            
            [psvEl,Workspace.Variables.aH2OEl] = electrolyteParameters(Workspace.Variables.T,p2,electrolyte);
            
            Workspace.Variables.ps = p1 - psvEl; % bar, Partial pressure of both hydrogen and oxygen in the system
            
            if any(Workspace.Variables.ps < 0)
                error('Pressure too low! System pressure lower than saturated vapor pressure of the electrolyte solution.')
            end
            
            % Nernst equation pressure correction
            funcHandle = @(Workspace) (Workspace.Constants.R.*Workspace.Variables.T)/(Workspace.Constants.n_e*Workspace.Constants.F).*log(Workspace.Variables.ps.^(3/2)./Workspace.Variables.aH2OEl);
    end
    
    Workspace.Coefficients = struct();
    
    Ucor = func(funcHandle,Workspace);
end