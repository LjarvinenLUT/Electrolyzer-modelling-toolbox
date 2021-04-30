% Activation overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article

function output = activation(varargin)
    
    defaultModel = 2;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
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
            Uact = @(coeffs,vars) 2*((R.*vars.T)./(n_e*F)).*asinh(vars.Current./(2*coeffs.j0));
            coeffs = struct('alpha',1/2,'j0',[]);
        case 2 % Hyperbolic sine approximation with variable alpha
            Uact = @(coeffs,vars) 1/coeffs.alpha.*((R.*vars.T)./(n_e*F)).*asinh(vars.Current./(2*coeffs.j0));
            coeffs = struct('alpha',[],'j0',[]);
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            Uact = @(coeffs,vars) 1/coeffs.alpha.*((R.*vars.T)./(n_e*F)).*log(vars.Current./coeffs.j0);
            coeffs = struct('alpha',[],'j0',[]);
    end

    output = struct('name','Uact','func',Uact,'coeffs',coeffs);
    
end