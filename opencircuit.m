% Open circuit voltage
% Inputs:   V - Measured values
%           type - Electrolysis type, "PEM" or "alkaline"
%           model - Used model reference, numeric, from which article

function Uocv = opencircuit(V,type,model)
    
    T = V(1);
    pH2 = V(2);
    pO2 = V(3);
    
    U_0 = reversible(T,model); %1.5184 -1.5421e-3.*T + 9.523e-5.*T.*log(T) + 9.84e-8.*T.^2;

    Uocv = U_0 + f*(pH2.*pO2.^(1/2));

end