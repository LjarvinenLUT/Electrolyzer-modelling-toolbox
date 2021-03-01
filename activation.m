% Activation overpotetial function handle generator
% Inputs:   V - Measured variables
%           D - Parameters
%           model - Used model reference, numeric, from which article

function Uact = activation(varargin)
    
    defaultModel = 1;

    parser = inputParser;
%     addRequired(parser,'T',@(x) isnumeric(x));
%     addRequired(parser,'j',@(x) isnumeric(x));
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    

    
    global F n_e R
        
    switch model
        case 1 % Hyperbolic sine approximation with alpha assumed to be 1/2
            Uact = @(j,T,j0) 2*((R*T)./(n_e*F)).*asinh(j./(2*j0));
        case 2 % Hyperbolic sine approximation with variable alpha
            Uact = @(j,T,j0,a) 1/a*((R*T)./(n_e*F)).*asinh(j./(2*j0));
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            Uact = @(j,T,j0,a) ((R*T)./(n_e*a*F)).*log(j./j0);
    end

end