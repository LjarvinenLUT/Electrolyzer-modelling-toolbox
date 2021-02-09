% Overpotentials

function U = overpotentials(T)
    
    Uocv = @opencircuit;
    Uact = @activation;
    Uohm = @ohmic;

    U = @(D) Uocv(T,D) + Uact(D) + Uohm(D);

end