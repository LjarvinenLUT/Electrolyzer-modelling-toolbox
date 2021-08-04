function [aH2OEl] = electrolyteWaterActivity(T,m,electrolyte,varargin)
% ELECTROLYTEWATERACTIVITY Calculate the water activity for a water 
%   solution of KOH or NaOH.
%
%   aH2OEl = ELECTROLYTEWATERACTIVITY(T,m,electrolyte) calculates values
%       using equations from Balej, J., "Water Vapour Partial Pressures and 
%       Water Activities in Potassium and Sodium Hydroxide Solutions Over 
%       Wide Concentration and Temperature Ranges", 
%       J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
%       - Inputs:   T -- Temperature in Kelvin
%                   m -- Solution molality, mol of solute/kg of solvent
%                   electrolyte -- Electrolyte type: 1 = KOH
%                                                    2 = NaOH
%       -Outputs:   aH2OEl -- Water activity for the solution
%
%   aH2OEl = ELECTROLYTEPARAMETERS(_,'model',model)
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

[~,aH2OEl] = electrolyteParameters(T,m,electrolyte,model);
