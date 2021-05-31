function Uocv = nernst(T,p1,p2,varargin)
% NERNST  Create a func object for calculation of open circuit potential
% using Nernst equation in modelling of water electrolysis.
%
%   Uocv = NERNST(T,p1,p2,'type',t) creates func object for
%           electrolyzer type defined by variable t. 
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
%   Uocv = NERNST(_,'type','alkaline','electrolyte',e) changes the 
%           electrolyte when alkaline electrolysis is considered.
%           Available electrolytes are:
%               KOH -- Potassium hydroxide (default)
%               NaOH -- Sodium hydroxide
%
%   Uocv = NERNST(_,'model',n) uses a model defined by number n for the
%           reversible voltage. Models available are presented in function 
%           reversible
% 
%   NERNST calls functions reversible and nernstPressureCorrection and
%   combines their output FUNC objects.
% 
%   See also ACTIVATION, OHMIC, CONCENTRATION, FUNC, REVERSIBLE,
%   NERNSTPRESSURECORRECTION

    %% Parse inputs
    defaultModel = 6;
    defaultType = 'NaN';
    defaultElectrolyte = 'KOH';

    parser = inputParser;
    addRequired(parser,'T',@(x) isnumeric(x));
    addRequired(parser,'p1',@(x) isnumeric(x));
    addRequired(parser,'p2',@(x) isnumeric(x));
    addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
    
    parse(parser,T,p1,p2,varargin{:});
    
    type = string(lower(parser.Results.type));
    if strcmp(type,"alkali")
        type = "alkaline";
    end
    model = parser.Results.model;
    electrolyte = string(parser.Results.electrolyte);
    
    %% Print parameters to command window
    fprintf('\nOpen circuit voltage calculation properties:\n')
    fprintf('Reversible voltage model: %d\n', model)
    fprintf('Electrolyzer type: %s\n', type)
    if strcmp(type,"alkaline")
        fprintf('Electrolyte: %s\n', electrolyte)
        if ~strcmp(electrolyte,"KOH")&&~strcmp(electrolyte,"NaOH")
            error("Only KOH and NaOH defined as alkaline electrolytes. Define the electrolyte by calling this function with parameter 'electrolyte' followed by one of the options.")
        end
    elseif ~strcmp(type,"pem")
        error("Only PEM and alkaline electrolysis defined for Nernst equation. Define the electrolyzer type by calling function with parameter 'type' followed by one of the options.")
    end
        
    
    
    %% Nernst equation
    
    Uocv = func.add(reversible(T,model),nernstPressureCorrection(T,p1,p2,'type',type,'electrolyte',electrolyte));
    
end