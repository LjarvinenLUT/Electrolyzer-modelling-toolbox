% Activation overpotetial
% Inputs:   V - Measured variables
%           D - Parameters
%           model - Used model reference, numeric, from which article

function Uact = activation(V, D, model)
    
    T = V(1);
    a = D(3);
    j0 = D(4);
    j = D(5);
    
    global F n_e R
        
    switch model
        case 1 % Buttler-Volmer equation
        case 2 % Buttler-Volmer equation, equal concentrations
        case 3 % Bockris et al. "Fuel cells: Their Electrochemistry", McGraw-Hill, New Your, 1969
            Uact = ((R*T)./F).*asinh(j./(2*j0));
        case 4 % Dale et al. "Semiempirical model based on thermodynamic principles for determining 6kW proton exchange membrane electrolyzer stack charactersitics"
            Uact = ((R*T)./(2*a_k*F)).*asinh(j./(2*j0));
        case 5 
            Uact = ((R*T)./(a_k*F)).*asinh(j./(2*j0));
        case 6 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            Uact = ((R*T)./(2*a_k*F)).*log(j./j0);
    end

end