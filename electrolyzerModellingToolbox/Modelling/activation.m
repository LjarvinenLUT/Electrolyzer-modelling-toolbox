function Uact = activation(varargin)
% ACTIVATION  Create a func object for activation overpotential in 
% modelling of water electrolysis.
%
%   Uact = ACTIVATION() uses hyperbolic sine approximation with variable
%   electron transfer coefficient, alpha
%
%   Uact = ACTIVATION('model',n) uses a model defined by number n.  
%
%   Models available currently are:
%       #1 -- Hyperbolic sine approximation with alpha assumed to be 1/2
%       #2 -- Hyperbolic sine approximation with variable alpha
%       #3 -- Tafel equation
%
%   All the models require the following func Workspace structure values:
%       Constants (obtained by the function automatically):
%           F -- Faraday's constant
%           R -- Universal gas constant
%           n_e -- Number of electrons transferred in a single
%                   electrochemical water splitting reaction
%       Variables: 
%           current -- Current density
%           T -- Temperature in kelvins
%   	Coefficients:
%           alpha -- Electron transfer coefficient
%           j0 -- Exchange current density
%
%   See also NERNST, OHMIC, CONCENTRATION, FUNC
    
    %% Parse input
    
    defaultModel = 2;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    
    %% Define the func object

    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = getConstants;
    
    Workspace.Variables = struct('current',[],'T',[]);
    
    switch model
        case 1 % Hyperbolic sine approximation with alpha assumed to be 1/2
            modelStr = model + " -- Hyperbolic sine approximation with alpha assumed to be 1/2";
            Workspace.Coefficients = struct('alpha',1/2,'j0',[]);
            Fitlims = struct('j0',{{1e-15,1e-5,1}});
            funcHandle = @(Workspace) 2*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*asinh(Workspace.Variables.current./(2*Workspace.Coefficients.j0));
        case 2 % Hyperbolic sine approximation with variable alpha
            modelStr = model + " -- Hyperbolic sine approximation with variable alpha";
            Workspace.Coefficients = struct('alpha',[],'j0',[]);
            Fitlims = struct('alpha',{{0,0.1,1}},'j0',{{1e-15,1e-5,1}});
            funcHandle = @(Workspace) 1/Workspace.Coefficients.alpha.*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*asinh(Workspace.Variables.current./(2*Workspace.Coefficients.j0));
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            modelStr = model + " -- Tafel equation";
            Workspace.Coefficients = struct('alpha',[],'j0',[]);
            Fitlims = struct('alpha',{{0,0.1,1}},'j0',{{1e-15,1e-5,1}});
            funcHandle = @(Workspace) 1/Workspace.Coefficients.alpha.*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log(Workspace.Variables.current./Workspace.Coefficients.j0);
        otherwise
            error("Activation overpotential model #" + num2str(model) + " not defined.")
    end
    

    Uact = func(funcHandle,Workspace,Fitlims);
    
    %% Print information to command window
    
    fprintf('\nActivation overpotential modelling properties:\n')
    % TODO: Switch statements could be maybe combined and printing done after the
    % functionality. However this prevents printing the settings if error
    % occurs
    fprintf('Model: %s\n', modelStr)

        
    
end