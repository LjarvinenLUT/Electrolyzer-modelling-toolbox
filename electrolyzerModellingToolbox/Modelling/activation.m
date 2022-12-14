function Uact = activation(varargin)
% ACTIVATION  Create a func object for activation overpotential in 
% modelling of water electrolysis.
%
%   Uact = ACTIVATION() uses hyperbolic sine approximation with variable
%   electron transfer parameter, alpha
%
%   Uact = ACTIVATION('model',n) uses a model defined by number n.  
%
%   Models available currently are:
%       #1 -- Hyperbolic sine approximation with variable alpha
%       #2 -- Hyperbolic sine approximation with alpha assumed to be 1/2
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
%           T -- Temperature in kelvin
%   	Parameters:
%           alpha -- Electron transfer coefficient
%           j0 -- Exchange current density
%
%   See also NERNST, OHMIC, CONCENTRATION, FUNC
    
    %% Parse input
    
    defaultModel = 1;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    
    %% Define the func object

    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = getConstants;
    
    Workspace.Variables = struct('current',[],'T',[]);
    
    Workspace.Dependencies.warning_j0 = "if ~isempty(Workspace.Parameters.j0)&&changed.T;"+...
                                "warning('Temperature dependency of the exchange current density (j0) is not taken into account in the models. Changing the temperature without a new fit will result in unrealistic behaviour.');"+...
                                "end;";
    
    switch model
        case 1 % Hyperbolic sine approximation with variable alpha
            modelStr = model + " -- Hyperbolic sine approximation with variable alpha";
            Workspace.Parameters = struct('alpha',[],'j0',[]);
            Fitlims = struct('alpha',{{0,0.1,1}},'j0',{{1e-15,1e-5,1}});
            funcHandle = @(Workspace) 1/Workspace.Parameters.alpha.*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*asinh(Workspace.Variables.current./(2*Workspace.Parameters.j0));
        case 2 % Hyperbolic sine approximation with alpha assumed to be 1/2
            modelStr = model + " -- Hyperbolic sine approximation with alpha assumed to be 1/2";
            Workspace.Parameters = struct('alpha',1/2,'j0',[]);
            Fitlims = struct('j0',{{1e-15,1e-5,1}});
            funcHandle = @(Workspace) 2*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*asinh(Workspace.Variables.current./(2*Workspace.Parameters.j0));
        case 3 % Tafel equation (valid when j/j0 > 4 https://doi.org/10.1016/j.jpowsour.2005.03.174)
            modelStr = model + " -- Tafel equation";
            Workspace.Parameters = struct('alpha',[],'j0',[]);
            Fitlims = struct('alpha',{{0,0.1,1}},'j0',{{1e-15,1e-5,1}});
            funcHandle = @(Workspace) 1/Workspace.Parameters.alpha.*((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log(Workspace.Variables.current./Workspace.Parameters.j0);
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