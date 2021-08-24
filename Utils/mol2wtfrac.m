function wtfrac = mol2wtfrac(m,M)
%MOL2WTFRAC Convert concentration value in 
%   molality (moles of solute/mass of solvent in kg) to
%   weight fraction (mass of solute/mass of solution).
%
%   Inputs:
%       m -- concentration as molality
%       M -- Molar mass of the electrolyte
% 
%    Outputs:
%       wtfrac -- concentration in weight fraction
% 
% See also WTFRAC2MOL

wtfrac = 1/(1+1/(m*M));

end