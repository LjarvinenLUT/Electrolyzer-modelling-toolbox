% Activation overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article

function Uact = activation(varargin)
    
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
    
    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = get_constants;
    
    Workspace.Variables = struct('current',[],'T',[]);
        
    switch model
        case 1 % Hyperbolic sine approximation with alpha assumed to be 1/2
            Workspace.Coefficients = struct('alpha',1/2,'j0',[]);
            funcHandle = @(Workspace) 2*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*asinh(Workspace.Variables.current./(2*Workspace.Coefficients.j0));
        case 2 % Hyperbolic sine approximation with variable alpha
            Workspace.Coefficients = struct('alpha',[],'j0',[]);
            funcHandle = @(Workspace) 1/Workspace.Coefficients.alpha.*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*asinh(Workspace.Variables.current./(2*Workspace.Coefficients.j0));
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            Workspace.Coefficients = struct('alpha',[],'j0',[]);
            funcHandle = @(Workspace) 1/Workspace.Coefficients.alpha.*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log(Workspace.Variables.current./Workspace.Coefficients.j0);
    end

    Uact = func(funcHandle,Workspace);
    
end