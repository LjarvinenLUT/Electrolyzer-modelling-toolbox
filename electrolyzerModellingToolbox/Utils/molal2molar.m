function Molarity = molal2molar(molality,T,solute,printWarning)
% MOLAL2MOLAR Convert concentration value in 
%   molality (moles of solute/mass of solvent in kg) to
%   molarity (moles of solute/volume of solution in L).
%
%   Inputs:
%       molality -- concentration in molality
%       T -- temperature in kelvin
%       solute -- chemical formula of the solute
%       printWarning -- true or false for either printing the warnings for
%           exceedin temperature and concentration limits of solution
%           density data or not.
% 
%    Outputs:
%       Molarity -- concentration as molarity
% 
% See also MOLAR2MOLAL, WTFRAC2MOLAL, MOLAL2WTFRAC

% If printWarning not determined
if ~exist("printWarning","var")
    printWarning = true;
end

% Molar mass of solute (g/mol)
switch solute
    case 'KOH'
        M = 39.0983 + 15.9994 + 1.0079;
    case 'NaOH'
        M = 22.9898 + 15.9994 + 1.0079;
end

Molarity = solutionDensity(T,molality,solute,printWarning)./( 1./molality + M*1e-3 );

end