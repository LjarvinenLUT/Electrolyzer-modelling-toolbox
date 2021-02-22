% Antoine equation for determining saturated water vapor pressure
% Inputs:   T - Measured temperature
%          

function psv = antoine(T)

if T<273||T>374
    error('Antoine equation defined for temperatures between 273 K and 374 K')
end

%% Parameters

Tl = [273 304;
    304 334;
    334 363;
    363 374]; % Temperature limits
A = [5.40221;
    5.20389;
    5.0768;
    5.08354];
B = [1836.675;
    1733.926;
    1659.793;
    1663.125];
C = [-31.737;
    -39.485;
    -45.854;
    -45.622];

i = T>=Tl(:,1)&T<Tl(:,2);

%% Equation

psv = 10^(A(i)-B(i)/(T+C(i)));

end