% Nernst equation for open circuit voltage
% Inputs:   T - Measured temperature
%           p1 - Parameter 1: 
%               for "PEM" cathode pressure, in bar
%               for "alkaline" system pressure, in bar
%           p2 - Parameter 2: 
%               for "PEM" anode pressure, in bar
%               for "alkaline" electrolyte molality, in mol/kg of solvent
%           type - Electrolysis type, "PEM" or "alkaline"
%           model - Used reversible potential model reference, numeric, from which article
%           electrolyte - Electrolyte used for "alkaline"

function Uocv = nernst(T,p1,p2,varargin)
    
    defaultModel = 6;
    defaultElectrolyte = 'KOH';

    parser = inputParser;
    addRequired(parser,'T',@(x) isnumeric(x));
    addRequired(parser,'p1',@(x) isnumeric(x));
    addRequired(parser,'p2',@(x) isnumeric(x));
    addParameter(parser,'type',@(x) ischar(x)||isstring(x));
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
    
    parse(parser,T,p1,p2,varargin{:});
    
    type = string(lower(parser.Results.type));
    model = parser.Results.model;
    electrolyte = string(parser.Results.electrolyte);
    
    %% Errors
    if strcmp(type,"alkaline")
        if ~strcmp(electrolyte,"KOH")&&~strcmp(electrolyte,"NaOH")
            error('Only KOH and NaOH defined as alkaline electrolytes')
        end
    elseif ~strcmp(type,"pem")
        error('Only PEM and alkaline electrolysis defined for Nernst equation.')
    end
    
    %% Nernst equation
    
    Uocv = reversible(T,model) + nerns_pressure_correction(T,p1,p2,type,electrolyte);
    
end