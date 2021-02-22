% Antoine equation for determining saturated water vapor pressure
% Inputs:   T - Measured temperature
%          

function psv = water_vapor_pressure(T,varargin)

defaultModel = 1; % Default experimental parameters

parser = inputParser;
addRequired(parser,'T',@(x) isnumeric(x));
addOptional(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))

parse(parser,T,varargin{:});

model = parser.Results.model;


%% Case structure

switch model
    case 1 % Antoine equation
        
        if any(T<273|T>374)
            error('Antoine equation defined for temperatures between 273 K and 374 K')
        end
        
        % Parameters, from NIST Chemistry WebBook SRD 69, https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185&Mask=4&Type=ANTOINE&Plot=on
        % Bridgeman, O.C.; Aldrich, E.W., Vapor Pressure Tables for Water, J. Heat Transfer, 1964, 86, 2, 279-286, https://doi.org/10.1115/1.3687121 
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
        
        psv = NaN(size(T));
        for j = 1:numel(T)
            i = T(j)>=Tl(:,1)&T(j)<Tl(:,2); % index
            
            % Equation
            psv(j) = 10^(A(i)-B(i)/(T(j)+C(i)));
        end
        
    case 2 % Bridgeman-Aldrich equation, from  Balej, J., "Water Vapour Partial Pressures and Water Activities in Potassium and Sodium Hydroxide Solutions Over Wide Concentration and Temperature Ranges", J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
        % Does not work!!!
        % Way too complicated... Adds just unnecessary complication
        % Parameters
        t = T - 273.15; % Temperature in degC
        A = 1.06994980;
        Y1 = 4.16385282*(t-187)./(t+237.098157);
        B = 1.0137921;
        C = 5.83531e-4;
        Z = -1.87 + 3.74*(1.152894 - 0.745794*acosh(654.2906./(t+266.778)));
        alpha = Z.^2.*(187^2 - Z.^2)./(0.30231574.*(1 + 3.377565e-3.*t));
        Y2 = 3^1.5/(2*1.87^3).*(0.01*(t-187-alpha)).*(187^2 - (0.01*(t-187-alpha)).^2);
        
        % Equation
        
        psv = 10.^(A + Y1 - B.*(1 + C.*t).*Y2);
        
    case 3 % Experimental fit: Balej, J., "Water Vapour Partial Pressures and Water Activities in Potassium and Sodium Hydroxide Solutions Over Wide Concentration and Temperature Ranges", J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
        
        psv = 10.^(35.4462 - 3343.93./T - 10.9.*log10(T) + 0.0041645.*T);
        
end

end