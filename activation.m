% Activation overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article

function Uact = activation(T,varargin)
    
    defaultModel = 2;

    parser = inputParser;
    addRequired(parser,'T',@(x) isnumeric(x))
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,T,varargin{:});
    
    model = parser.Results.model;
    
    % Print parameters to command window
    fprintf('\nActivation overpotential calculation properties:\n')
    % TODO: Switch statements could be maybe combined and printing done after the
    % functionality. However this prevents printing the settings if error
    % occurs
    switch model
        case 1
            model_str = "Hyperbolic sine approximation with alpha assumed to be 1/2";
        case 2
            model_str = "Hyperbolic sine approximation with variable alpha";
        case 3
            model_str = "Tafel equation";
    end
    fprintf('Activation voltage model: %s\n', model_str)
    
    [F,R,n_e] = get_constants;
        
    switch model
        case 1 % Hyperbolic sine approximation with alpha assumed to be 1/2
            Uact = @(j0,a,j) 2*((R.*T)./(n_e*F)).*asinh(j./(2*j0));
        case 2 % Hyperbolic sine approximation with variable alpha
            Uact = @(j0,a,j) 1/a.*((R.*T)./(n_e*F)).*asinh(j./(2*j0));
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            Uact = @(j0,a,j) 1/a.*((R.*T)./(n_e*F)).*log(j./j0);
    end

end