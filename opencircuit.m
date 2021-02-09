% Open circuit voltage

function Uocv = opencircuit(T,D)
    
%     T = D(1);
    pH2 = D(1);
    pO2 = D(2);
    
    U_0 = 1.5184 -1.5421e-3.*T + 9.523e-5.*T.*log(T) + 9.84e-8.*T.^2;

    Uocv = U_0 + (pH2.*pO2.^(1/2));

end