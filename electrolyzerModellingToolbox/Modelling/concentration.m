function Ucon = concentration(varargin)
% CONCENTRATION  Create a func object for concentration overpotential in 
% modelling of water electrolysis.
%
%   Ucon = CONCENTRATION() uses concentration overpotential model with
%   logarithmic dependency from the fraction between current and its
%   limiting value, j_lim (model #2).
%
%   Ucon = CONCENTRATION('model',n) uses a model defined by number n.  
%
%   Models available currently are:
%       #1 -- Model with logarithmic dependency from the fraction between 
%               current and its limiting value.
%               Parameter:
%                   j_lim -- Limiting current density
%               Variables:
%                   current -- Current density
%                   T -- Temperature in kelvins
%       #2 -- Model with dependency on reagent concentrations on both anode
%               and cathode before and after mass transfer. This model
%               cannot be applied unless the reagent concentrations have
%               been measured. The model has no fitting parameters.
%               Recured variables are:
%                   T -- Temperature in kelvins
%                   Cx10 -- Equilibrium reagent concentration at the cathode
%                   Cx11 -- Reagent concentration at the cathode after mass
%                           transfer
%                   Cx20 -- Equilibrium reagent concentration at the anode
%                   Cx21 -- Reagent concentration at the anode after mass
%                           transfer
%
%   All the models require the following func Workspace structure values:
%       Constants (obtained by the function automatically):
%           F -- Faraday's constant
%           R -- Universal gas constant
%           n_e -- Number of electrons transferred in a single
%                   electrochemical water splitting reaction
%
%   See also NERNST, ACTIVATION, OHMIC, FUNC
    
%% Parse input
    defaultModel = 1;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    
    
    %% Constants
    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = getConstants;
    
    fprintf('\nConcentration overpotential modelling properties:\n')
    

    %% Define the func object
    
    switch model
        case 1 % Model with logarithmic dependency from the fraction between current and its limiting value.
            modelStr = model + " -- Limiting current model";
            funcHandle = @(Workspace) -((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log(1 - Workspace.Variables.current/Workspace.Parameters.j_lim);
            Workspace.Parameters = struct('j_lim',[]);
            Workspace.Variables = struct('current',[],'T',[]);
            Fitlims.j_lim = {'max(x)','max(x)*1.01','max(x)*3'};
        case 2 % F. Marangio "Theoretical model and experimental analysis of a high pressure PEM water electrolyzer for hydrogen production" https://doi.org/10.1016/j.ijhydene.2008.11.083
            % Cxn0 - reagent X concentration at the membrane-electrode
            %   interface, n = 1 for cathode, 2 for anode
            % Cxn1 - Concentration after mass transfer has occurred, 
            %   n = 1 for cathode, 2 for anode
            modelStr = model + " -- Model based on reagent concentrations on electrode surfaces";
            funcHandle = @(Workspace) ((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log((Workspace.Variables.Cx11/Workspace.Variables.Cx10)*(Workspace.Variables.Cx21/Workspace.Variables.Cx20)^(1/2)); % Both concentrations are functions of current
            Workspace.Variables = struct('T',[],'Cx11',[],'Cx10',[],'Cx21',[],'Cx20',[]);
            Workspace.Parameters = struct();
            Fitlims = struct([]);
            warning("Concentration overpotential model #" + num2str(model) + " cannot be used for fitting! Concentration on electrode surface is an unknown function of current. The function handle combines effects of hydrogen and oxygen side. Refer to manual for more information")
        otherwise
            error("Concentration overpotential model #" + num2str(model) + " not defined.")
    end
    fprintf('Model: %s\n', modelStr)
    
    Ucon = func(funcHandle,Workspace,Fitlims);
    
end