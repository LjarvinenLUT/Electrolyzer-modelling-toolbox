function Ucor = nernstPressureCorrection(type)
% NERNSTPRESSURECORRECTION  Create a func object for pressure correction 
% term of Nernst equation in modelling of water electrolysis.
%
%   Ucor = NERNSTPRESSURECORRECTION(type) creates func object
%           for electrolyzer type defined. 
%           Electrolyzer types availabel are:
%               PEM -- Polymer electrolyte membrane electrolysis
%               Alkaline -- Alkaline electrolysis.
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
%           psv -- Saturated water vapor pressure (for PEM only; dependency
%               on temperature)
%           ps -- System pressure (for alkaline only)
%           m -- Electrolyte molality (for alkaline only)
%           aH2OEl -- Water activity in electrolyte solution (for alkaline
%               only; dependency on temperature, molality and the
%               electrolyte)
%           psvEl -- Saturated water vapor pressure for electrolyte
%               solution (for alkaline only; dependency on temperature,
%               molality and the electrolyte)
% 
%   See also NERNST, REVERSIBLE, FUNC, WATERVAPORPRESSURE,
%   ELECTROLYTEPARAMETERS
    
%% Parse input

    type = string(lower(type));
    if strcmp(type,"alkali")
        type = "alkaline";
    end
    
    %% Errors
    if ~strcmp(type,"alkaline")&&~strcmp(type,"pem")
        error('Only PEM and alkaline electrolysis defined for Nernst equation.')
    end
    %% Constants
    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = getConstants;

    %% Define the func object

    Workspace.Variables.T = [];
    
    switch type
        case "pem"
            % Vapor pressure from Antoine equation
            Workspace.Dependencies.psv = 'Workspace.Variables.psv = waterVaporPressure(Workspace.Variables.T);'; % Dependency format
%             eval(Workspace.Dependencies.psv);
            
            Workspace.Variables.pCat = []; % bara, Cathode pressure
            
            Workspace.Variables.pAn = []; % bara, Anode pressure
            
%             if any(Workspace.Variables.pCat < Workspace.Variables.psv | Workspace.Variables.pAn < Workspace.Variables.psv)
%                 error('Pressure too low! Anode or cathode pressure lower than the saturated vapor pressure at the given temperature. ')
%             end
                      
            % Nernst equation pressure correction
            funcHandle = @(Workspace) (Workspace.Constants.R.*Workspace.Variables.T)/(Workspace.Constants.n_e*Workspace.Constants.F).*log((Workspace.Variables.pCat-Workspace.Variables.psv).*(Workspace.Variables.pAn-Workspace.Variables.psv).^(1/2));
            
        case "alkaline"
            
            
            Workspace.Variables.ps = []; % bara, System pressure
            
            Workspace.Variables.m = []; % mol/kg of solvent, Electrolyte molality

            Workspace.Variables.electrolyte = []; % 1 = KOH, 2 = NaOH, helper variable for determining the electrolyte related parameters
            
            % Electrolyte parameters
            Workspace.Dependencies.electrolyteParameters = '[Workspace.Variables.psvEl,Workspace.Variables.aH2OEl] = electrolyteParameters(Workspace.Variables.T,Workspace.Variables.m,Workspace.Variables.electrolyte);'; % Dependency format
%             eval(Workspace.Dependencies.electrolyteParameters);
            
%             if any(Workspace.Variables.ps < Workspace.Variables.psvEl)
%                 error('Pressure too low! System pressure lower than saturated vapor pressure of the electrolyte solution.')
%             end
            
            % Nernst equation pressure correction
            funcHandle = @(Workspace) (Workspace.Constants.R.*Workspace.Variables.T)/(Workspace.Constants.n_e*Workspace.Constants.F).*log((Workspace.Variables.ps-Workspace.Variables.psvEl).^(3/2)./Workspace.Variables.aH2OEl);
    end
    
    Workspace.Coefficients = struct();
    
    Ucor = func(funcHandle,Workspace);
end