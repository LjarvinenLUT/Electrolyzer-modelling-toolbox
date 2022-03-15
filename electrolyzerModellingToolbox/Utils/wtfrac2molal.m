function molality = wtfrac2molal(wtfrac,M)
% WTFRAC2MOLAL Convert concentration value in 
%   weight fraction (mass of solute/mass of solution) to 
%   molality (moles of solute/mass of solvent in kg).
%
%   Inputs:
%       wtfrac -- concentration in weight fraction, or percentage which is
%           assumed to be used if wtfrac is greater than 1. If wtfrac is 
%           given as percentage, it is automatically converted to weight 
%           fraction.
%       M -- Molar mass of the solute (g/mol)
% 
%    Outputs:
%       molality -- concentration as molality
% 
% See also MOLAL2WTFRAC, MOLAL2MOLAR, MOLAR2MOLAL


if wtfrac > 1 % If weight frqaction is assumed to be given as percents
    wtfrac = wtfrac/100;
end

molality = 1./(M*1e-3*(1./wtfrac-1));

end