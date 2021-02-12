% Open circuit voltage
% Inputs:   V - Measured values
%           type - Electrolysis type, "PEM" or "alkaline"
%           model - Used model reference, numeric, from which article

function Uocv = opencircuit(V,type,model)
    
    global F n_e R;

    T = V(1);
    pH2 = V(2);
    pO2 = V(3);
    
    % Reversible voltage
    U_0 = reversible(T,model);
    
    % Nernst equation
    Uocv = U_0 + (R.*T)/(n_e*F)*log(pH2.*pO2.^(1/2));

end