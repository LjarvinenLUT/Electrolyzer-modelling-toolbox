% Activation overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article

function Uact = activation(varargin)
    
    defaultModel = 2;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    

    
    [F,R,n_e] = get_constants;
        
    switch model
        case 1 % Hyperbolic sine approximation with alpha assumed to be 1/2
            Uact = @(j0,a,T,j) 2*((R.*T)./(n_e*F)).*asinh(j./(2*j0));
        case 2 % Hyperbolic sine approximation with variable alpha
            Uact = @(j0,a,T,j) 1/a.*((R.*T)./(n_e*F)).*asinh(j./(2*j0));
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            Uact = @(j0,a,T,j) 1/a.*((R.*T)./(n_e*F)).*log(j./j0);
    end

end