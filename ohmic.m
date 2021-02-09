% Ohmic overpotetial

function Uohm = ohmic(D)
    
    r = D(6);
    j = D(5);
    
    Uohm = r.*j;

end