% Activation overpotetial

function Uact = activation(D)
    
    a = D(3);
    j0 = D(4);
    j = D(5);
    
    Uact = a.*log(j./j0);

end