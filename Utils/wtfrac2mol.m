function m = wtfrac2mol(wtfrac,M)
% WTFRAC2MOL Convert concentration value in 
%   weight fraction (mass of solute/mass of solution) to 
%   molality (moles of solute/mass of solvent in kg).
%
%   Inputs:
%       wtfrac -- concentration in weight fraction, or percentage which is
%           assumed to be used if wtfrac is greater than 1. If wtfrac is 
%           given as percentage, it is automatically converted to weight 
%           fraction.
%       M -- Molar mass of the electrolyte
% 
%    Outputs:
%       m -- concentration as molality
% 
% See also MOL2WTFRAC


if wtfrac > 1 % If weight frqaction is assumed to be given as percents
    wtfrac = wtfrac/100;
end

m = 1/(M*(1/wtfrac-1));

end