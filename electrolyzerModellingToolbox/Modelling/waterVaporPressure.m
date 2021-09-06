function psv = waterVaporPressure(T,varargin)
% WATERVAPORPRESSURE Calculate water vapor pressure in temperature T in
% kelvins. Output is in units of bars.
%
%   psv = WATERVAPORPRESSURE(T) calculates water vapor pressure using
%       equation presented by Balej, J. in "Water Vapour Partial Pressures 
%       and Water Activities in Potassium and Sodium Hydroxide Solutions 
%       Over Wide Concentration and Temperature Ranges", 
%       J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
%
%   psv = WATERVAPORPRESSURE(T,'model',m) uses one of the two available
%       models:
%           #1 -- Antoine equation [bara]
%           #2 -- Equation presented by Balej, J. (default) [bara]
%
% See also ELECTROLYTEPARAMETERS

defaultModel = 2; % Default experimental parameters

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
        
    case 2 % Experimental fit: Balej, J., "Water Vapour Partial Pressures and Water Activities in Potassium and Sodium Hydroxide Solutions Over Wide Concentration and Temperature Ranges", J. Hydrogen Energy, Vol. 10, No. 4. pp. 233--243, 1985
        
        psv = 10.^(35.4462 - 3343.93./T - 10.9.*log10(T) + 4.1645e-3.*T);

end

end