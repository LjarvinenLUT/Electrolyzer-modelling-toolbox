% Concentration overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article
%               p_O2 - Oxygen partial pressure

function Ucon = concentration(varargin)
    
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
    [Workspace.Constants.F,Workspace.Constants.R,Workspace.Constants.n_e] = get_constants;
    
    fprintf('\nConcentration overpotential calculation properties:\n')
    fprintf('Electrolyzer: %d\n', model)
    %% Error checking
    
    if model == 3 && isempty(T)
        error('Temperature has to be given as a parameter for concentration overpotential function when model #3 is used.')
    end
    if model == 3 && isempty(p_O2)
        error('Variable "p_O2" (oxigen partial pressure at the anode) has to be given as a parameter for concentration overpotential function when model #3 is used.')
    end
    
    %%
    
    % Cx0 - reagent X concentration at the membrane-electrode interface
    % Cx1 - Concentration after mass transfer has occurred 
    % p_O2 - oxygen partial pressure (bar)
    % p_sat - vapor saturation pressure (bar)
    % T - Temperature (K)
    % b - Empiric coefficient
    % b2 - Concentration overpotential constant
    % j_lim - Current limiting term
    % j_max - Maximum current density (2 A/cm usually used)
    switch model
        case 1 % F. Marangio "Theoretical model and experimental analysis of a high pressure PEM water electrolyzer for hydrogen production" https://doi.org/10.1016/j.ijhydene.2008.11.083
            funcHandle = @(Workspace) ((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log((Workspace.Variables.Cx11/Workspace.Variables.Cx10)*(Workspace.Variables.Cx21/Workspace.Variables.Cx20)^(1/2)); % Both concentrations are functions of current
            Workspace.Variables = struct('T',[],'Cx11',[],'Cx10',[],'Cx21',[],'Cx20',[]);
            warning('Chosen concentration overpotential method cannot be used for fitting! Concentration on electrode surface is an unknown function of current. The function handle combines effects of hydrogen and oxygen side. Refer to manual for more information')
        case 2 % ???
            funcHandle = @(Workspace) -((Workspace.Constants.R.*Workspace.Variables.T)./(Workspace.Constants.n_e*Workspace.Constants.F)).*log(1 - Workspace.Variables.current/Workspace.Coefficients.j_lim);
            Workspace.Coefficients = struct('j_lim',[]);
            Workspace.Variables = struct('current',[],'T',[]);
        case 3 % J. Pukrushpan "Modeling and control for PEM fuel cell stack system" https://doi.org/10.1109/ACC.2002.1025268
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
            warning('Concentration overpotential model is based on fuel cell research where the amount of oxygen supplied to the cathode is limiting. Includes experimental numeric parameters that probably cannot be extended for electrolyzers!')
    end
    
    Ucon = func(funcHandle,Workspace);
    
end