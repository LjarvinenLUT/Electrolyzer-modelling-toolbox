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
%       #1 -- Model with dependency on reagent concentrations on both anode
%               and cathode before and after mass transfer. This model
%               cannot be applied unless the reagent concentrations have
%               been measured. The model has no fitting coefficients.
%               Recured variables are:
%                   T -- Temperature in kelvins
%                   Cx10 -- Equilibrium reagent concentration at the cathode
%                   Cx11 -- Reagent concentration at the cathode after mass
%                           transfer
%                   Cx20 -- Equilibrium reagent concentration at the anode
%                   Cx21 -- Reagent concentration at the anode after mass
%                           transfer              
%       #2 -- Model with logarithmic dependency from the fraction between 
%               current and its limiting value.
%               Coefficient:
%                   j_lim -- Limiting current density
%               Variables:
%                   current -- Current density
%                   T -- Temperature in kelvins
%       #3 -- Model obtained from fuel cell research that includes multiple
%               experimental parameters. Applicability of the model is
%               questionable.
%               Requires following parameter inputs as name-value pairs:
%                   T -- Temperature in kelvins
%                   p_O2 -- Oxygen partial pressure in bars
%               Coefficients:
%                   j_max -- Maximum current density
%                   b2 -- Concentration overpotential constant
%               Variables:
%                   current -- Current density
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
    defaultModel = 2;
    defaultTemperature = [];
    defaultP_O2 = [];

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    addParameter(parser,'Temperature',defaultTemperature,@(x) isnumeric(x));
    addParameter(parser,'p_O2',defaultP_O2,@(x) isnumeric(x));
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    T = parser.Results.Temperature; % This value is not directly inserted to Workspace.Variables because the value is only used in model 3 before the function handle.
    p_O2 = parser.Results.p_O2;
    
    
    %% Constants
    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = getConstants;
    
    fprintf('\nConcentration overpotential calculation properties:\n')
    fprintf('Electrolyzer: %d\n', model)
    %% Error checking
    
    if model == 3 && isempty(T)
        error('Temperature has to be given as a parameter for concentration overpotential function when model #3 is used.')
    end
    if model == 3 && isempty(p_O2)
        error('Variable "p_O2" (oxigen partial pressure at the anode) has to be given as a parameter for concentration overpotential function when model #3 is used.')
    end
    
    %% Define the func object
    
    switch model
        case 1 % F. Marangio "Theoretical model and experimental analysis of a high pressure PEM water electrolyzer for hydrogen production" https://doi.org/10.1016/j.ijhydene.2008.11.083
            % Cxn0 - reagent X concentration at the membrane-electrode
            %   interface, n = 1 for cathode, 2 for anode
            % Cxn1 - Concentration after mass transfer has occurred, 
            %   n = 1 for cathode, 2 for anode
            funcHandle = @(Workspace) ((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log((Workspace.Variables.Cx11/Workspace.Variables.Cx10)*(Workspace.Variables.Cx21/Workspace.Variables.Cx20)^(1/2)); % Both concentrations are functions of current
            Workspace.Variables = struct('T',[],'Cx11',[],'Cx10',[],'Cx21',[],'Cx20',[]);
            Workspace.Coefficients = struct();
            warning("Concentration overpotential model #" + num2str(model) + " cannot be used for fitting! Concentration on electrode surface is an unknown function of current. The function handle combines effects of hydrogen and oxygen side. Refer to manual for more information")
        case 2 % ???
            funcHandle = @(Workspace) -((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log(1 - Workspace.Variables.current/Workspace.Coefficients.j_lim);
            Workspace.Coefficients = struct('j_lim',[]);
            Workspace.Variables = struct('current',[],'T',[]);
        case 3 % J. Pukrushpan "Modeling and control for PEM fuel cell stack system" https://doi.org/10.1109/ACC.2002.1025268
            % p_O2 - oxygen partial pressure (bar)
            % p_sat - vapor saturation pressure (bar)
            % b - Empiric coefficient
            p_sat = water_vapor_pressure(T);
            px = p_O2/0.1173 + p_sat;
            atm = 1.01325;
            if px < 2 * atm
                b1 = (7.16*10^-4*T - 0.622)*px + (-1.45*10^-3*T + 1.68);
            else
                b1 = (8.66*10^-5*T - 0.068)*px + (-1.6*10^-4*T + 0.54);
            end
            
            funcHandle = @(Workspace) Workspace.Variables.current*(b1*Workspace.Variables.current/Workspace.Coefficients.j_max)^Workspace.Coefficients.b2;
            Workspace.Coefficients = struct('j_max',[],'b2',[]);
            Workspace.Variables = struct('T',T,'current',[]);
            
            warning("Concentration overpotential model #" + num2str(model) + " is based on fuel cell research where the amount of oxygen supplied to the cathode is limiting. Includes experimental numeric parameters that probably cannot be extended for electrolyzers!")
        otherwise
            error("Concentration overpotential model #" + num2str(model) + " not defined.")
    end
    
    Ucon = func(funcHandle,Workspace);
    
end