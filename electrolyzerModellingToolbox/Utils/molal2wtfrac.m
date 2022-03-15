function wtfrac = molal2wtfrac(molality,M)
% MOLAL2WTFRAC Convert concentration value in 
%   molality (moles of solute/mass of solvent in kg) to
%   weight fraction (mass of solute/mass of solution).
%
%   Inputs:
%       molality -- concentration as molality
%       M -- Molar mass of the solute (g/mol)
% 
%    Outputs:
%       wtfrac -- concentration in weight fraction
% 
% See also WTFRAC2MOLAL, MOLAL2MOLAR, MOLAR2MOLAL

wtfrac = 1./(1+1./(molality*M*1e-3));

end