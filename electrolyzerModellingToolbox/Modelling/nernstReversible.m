function Urev = nernstReversible(model)
% NERNSTREVERSIBLE  Create a func object for calculation of reversible potential
% of the water splitting reaction for Nernst equation.
% 
%   Urev = NERNSTREVERSIBLE(model) uses a model defined by parameter model,
%       provided as an integer
% 
%   NERNSTREVERSIBLE has six (6) models available. Their more in detail
%   description can be found in the mathematical documentation of this
%   tool.
%
%   The output func object has the following fields in workspace
%   structures:
%       Constants (obtained by the function automatically):
%           F -- Faraday's constant
%           R -- Universal gas constant
%           n_e -- Number of electrons transferred in a single
%                   electrochemical water splitting reaction
%       Variables:
%           T -- Temperature in kelvin
% 
%   See also NERNST, FUNC, NERNSTPRESSURECORRECTION,


[Workspace.Constants.F,~,Workspace.Constants.n_e] = getConstants;

Workspace.Variables = struct('T',[]);
    
    switch model
        case 1 % Koponen et al. "Effect of Power Quality on the Design of PEM Water Electrolysis Systems", Applied Energy 279, 115791, 2020, https://doi.org/10.1016/j.apenergy.2020.115791
            funcHandle = @(Workspace) 1/(Workspace.Constants.n_e*Workspace.Constants.F)*(-159.6.*Workspace.Variables.T + 2.8472e5);
        case 2 % Le Roy et al. "The thermodynamics of aqueous water electrolysis." Journal of the Electrochemical Society, 127:1954–1962, 1980. http://dx.doi.org/10.1149/1.2130044
            funcHandle = @(Workspace) 1.5184 - 1.5421e-3.*Workspace.Variables.T + 9.523e-5.*Workspace.Variables.T.*log(Workspace.Variables.T) + 9.84e-8.*Workspace.Variables.T.^2;
        case 3 % Dale et al. "Semiempirical model based on thermodynamic principles for determining 6 kW proton exchange membrane electrolyzer stack characteristics", Journal of Power Sources 185, Issue 2, 1348--1353, 2008, https://doi.org/10.1016/j.jpowsour.2008.08.054
            % Same as Le Roy, but with own measurements and curve fit...
            % Very inaccurate at high temperatures (error 5 mV)
            funcHandle = @(Workspace) 1.5241 - 1.2261e-3.*Workspace.Variables.T + 1.1858e-5.*Workspace.Variables.T.*log(Workspace.Variables.T) + 5.6692e-7.*Workspace.Variables.T.^2;
        case 4 % Hammoudi et al., "New multi-physics approach for modelling and design of alkaline electrolyzers", Int. Journ. of Hydr. Ene. 37, Issue 19, 13895--13913, 2012, https://doi.org/10.1016/j.ijhydene.2012.07.015
            funcHandle = @(Workspace) 1.50342 - 9.956e-4.*Workspace.Variables.T + 2.5e-7.*Workspace.Variables.T.^2;
        case 5 % D.M. Bernardi, M.W. Verbrugge, "A Mathematical Model of the Solid-Polymer-Electrolyte Fuel Cell", J. Electrochem. Soc., 139 (1992), p. 2477
            % Inaccurate at higher temperatures (error 1 mV)
            funcHandle = @(Workspace) 1.229 - 8.46e-4.*(Workspace.Variables.T - 298.15);
        case 6 % da Costa Lopes et al. "Experimental and theoretical development of a PEM electrolyzer model applied to energy storage systems", 2009 Brazilian Power Electronics Conference, DOI: 10.1109/COBEP.2009.5347619
            % Very inaccurate at higher temperatures (error 3.5 mV)
            funcHandle = @(Workspace) 1.449 - 6.139e-4.*Workspace.Variables.T - 4.592e-7.*Workspace.Variables.T.^2 + 1.46e-10.*Workspace.Variables.T.^3;
        otherwise
            error('No valid model number given for reversible voltage.')
    end
    
    Workspace.Coefficients = struct();
    
    Urev = func(funcHandle,Workspace);

end