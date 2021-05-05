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

function output = nernst(varargin)
    
    defaultModel = 6;
    defaultType = 'NaN';
    defaultElectrolyte = 'KOH';

    parser = inputParser;
    addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
    
    parse(parser,varargin{:});
    
    type = string(lower(parser.Results.type));
    model = parser.Results.model;
    electrolyte = string(parser.Results.electrolyte);
    
    % Print parameters to command window
    fprintf('\nOpen circuit voltage calculation properties:\n')
    fprintf('Electrolyzer type: %s\n', type)
    if strcmp(type,"alkaline")
        fprintf('Electrolyte: %s\n', electrolyte)
    end
    fprintf('Reversible voltage model: %d\n', model)
    
    %% Errors
    if strcmp(type,"alkaline")
        if ~strcmp(electrolyte,"KOH")&&~strcmp(electrolyte,"NaOH")
            error("Only KOH and NaOH defined as alkaline electrolytes. Define the electrolyte by calling this function with parameter 'electrolyte' followed by one of the options.")
        end
    elseif ~strcmp(type,"pem")
        error("Only PEM and alkaline electrolysis defined for Nernst equation. Define the electrolyzer type by calling function with parameter 'type' followed by one of the options.")
    end
    
    %% Nernst equation
    
    Uocv = @(coeffs,vars) reversible(vars.T,model) + nernst_pressure_correction(vars,'type',type,'electrolyte',electrolyte);
    coeffs = struct([]);
    
    output = struct('name','Uocv','func',Uocv,'coeffs',coeffs);
end