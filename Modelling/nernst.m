function Uocv = nernst(type,varargin)
% NERNST  Create a func object for calculation of open circuit potential
% using Nernst equation in modelling of water electrolysis.
%
%   Uocv = NERNST(type) creates func object for
%           electrolyzer type defined by variable type. 
%           Electrolyzer types availabel are:
%               PEM -- Polymer electrolyte membrane electrolysis
%               Alkaline -- Alkaline electrolysis.
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

    parser = inputParser;
    addRequired(parser,'type',@(x) ischar(x)||isstring(x));
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,type,varargin{:});
    
    type = string(lower(type));
    if strcmp(type,"alkali")
        type = "alkaline";
    end
    model = parser.Results.model;
    
    %% Print parameters to command window
    fprintf('\nOpen circuit voltage calculation properties:\n')
    fprintf('Reversible voltage model: %d\n', model)
    fprintf('Electrolyzer type: %s\n', type)
    if ~strcmp(type,"alkaline") && ~strcmp(type,"pem")
        error("Only PEM and alkaline electrolysis defined for Nernst equation. Define the electrolyzer type by inputing one of the options as the first parameter to nernst.")
    end
        
    
    
    %% Nernst equation
    
    Uocv = func.add(reversible(model),nernstPressureCorrection(type));
    
end