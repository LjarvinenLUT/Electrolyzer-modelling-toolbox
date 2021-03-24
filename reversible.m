% Reversible voltage
% Inputs:   T - Measured temperature
%           model - Used model reference, numeric, from which article

function U_0 = reversible(T,model)

[F,~,n_e] = get_constants;
    
    switch model
        case 1 % Le Roy et al. "The thermodynamics of aqueous water electrolysis." Journal of the Electrochemical Society, 127:1954–1962, 1980. http://dx.doi.org/10.1149/1.2130044
            U_0 = 1.5184 - 1.5421e-3.*T + 9.523e-5.*T.*log(T) + 9.84e-8.*T.^2;
        case 2 % Dale et al. "Semiempirical model based on thermodynamic principles for determining 6 kW proton exchange membrane electrolyzer stack characteristics", Journal of Power Sources 185, Issue 2, 1348--1353, 2008, https://doi.org/10.1016/j.jpowsour.2008.08.054
            % Same as Le Roy, but with own measurements and curve fit...
            % Very inaccurate at high temperatures (5 mV)
            U_0 = 1.5241 - 1.2261e-3.*T + 1.1858e-5.*T.*log(T) + 5.6692e-7.*T.^2;
        case 3 % Hammoudi et al., "New multi-physics approach for modelling and design of alkaline electrolyzers", Int. Journ. of Hydr. Ene. 37, Issue 19, 13895--13913, 2012, https://doi.org/10.1016/j.ijhydene.2012.07.015
            U_0 = 1.50342 - 9.956e-4.*T + 2.5e-7.*T.^2;
        case 4 % D.M. Bernardi, M.W. Verbrugge, "A Mathematical Model of the Solid-Polymer-Electrolyte Fuel Cell", J. Electrochem. Soc., 139 (1992), p. 2477
            % Inaccurate at higher temperatures (1 mV)
            U_0 = 1.229 - 8.46e-4.*(T - 298.15);
        case 5 % Costa Lopes et al. "Experimental and theoretical development of a PEM electrolyzer model applied to energy storage systems", 2009 Brazilian Power Electronics Conference, DOI: 10.1109/COBEP.2009.5347619
            % Very inaccurate at higher temperatures (3.5 mV)
            U_0 = 1.449 - 6.139e-4.*T - 4.592e-7.*T.^2 + 1.46e-10.*T.^3;
        case 6 % Koponen et al. "Effect of Power Quality on the Design of PEM Water Electrolysis Systems", Applied Energy 279, 115791, 2020, https://doi.org/10.1016/j.apenergy.2020.115791
            U_0 = 1/(n_e*F)*(-159.6.*T + 2.8472e5);
    end

end