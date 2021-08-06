function [psvEl] = electrolyteWaterVaporPressure(T,m,electrolyte,varargin)
% ELECTROLYTEWATERVAPORPRESSURE Calculate saturated water vapor pressure 
%   for a water solution of KOH or NaOH.
%
%   psvEl = ELECTROLYTEWATERVAPORPRESSURE(T,m,electrolyte) calculates values
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
%
%   psvEl = ELECTROLYTEWATERVAPORPRESSUR(_,'model',model)
%       allows changing the model for one of the two contained:
%           #1 -- Equations presented by Balej J. (default)
%           #2 -- Ideal solution approximation
%
% See also ELECTROLYTEPARAMETERS

defaultModel = 1; % Default experimental parameters

parser = inputParser;
addRequired(parser,'T',@(x) isnumeric(x)); % K, Temperature
addRequired(parser,'m',@(x) isnumeric(x)); % mol/kg of solvent, Molality
addRequired(parser,'electrolyte',@(x) isnumeric(x)&&isscalar(x)); % Electrolyte, 1 = KOH, 2 = NaOH
addOptional(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))

parse(parser,T,m,electrolyte,varargin{:});

model = parser.Results.model;

[psvEl,~] = electrolyteParameters(T,m,electrolyte,model);
