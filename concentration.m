% Concentration overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article
%               p_O2 - Oxygen partial pressure

function output = concentration(varargin)
    
    defaultModel = 2;

    parser = inputParser;
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    model = parser.Results.model;
    
    %% Global parameters
    [F,R,n_e] = get_constants();

    f = R/(n_e*F);
    
    fprintf('\nConcentration overpotential calculation properties:\n')
    fprintf('Electrolyzer: %d\n', model)
    %% Error checking
    
    if model == 3 && ~isnumeric(vars.p_O2)
        error('Variable "p_O2" (oxigen partial pressure at the anode) has to be set for concentration model 3')
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
            Ucon = @(vars) f*vars.T*log((vars.Cx11/vars.Cx10)*(vars.Cx21/vars.Cx20)^(1/2)); % Both concentrations are functions of current
            coeffs = struct();
            warning('Chosen concentration overpotential method cannot be used for fitting! Concentration on electrode surface is an unknown function of current. The function handle combines effects of hydrogen and oxygen side. Refer to manual for more information')
        case 2 % ???
            Ucon = @(coeffs,vars) -f*vars.T*log(1 - vars.Current/coeffs.j_lim);
            coeffs = struct('j_lim',[]);
        case 3 % J. Pukrushpan "Modeling and control for PEM fuel cell stack system" https://doi.org/10.1109/ACC.2002.1025268
            p_sat = water_vapor_pressure(vars.T);
            px = vars.p_O2/0.1173 + p_sat;
            atm = 1.01325;
            if px < 2 * atm
                b1 = (7.16*10^-4*vars.T - 0.622)*px + (-1.45*10^-3*vars.T + 1.68); 
            else
                b1 = (8.66*10^-5*vars.T - 0.068)*px + (-1.6*10^-4*vars.T + 0.54);
            end
            
            Ucon = @(coeffs,vars) vars.Current*(b1*vars.Current/coeffs.j_max)^coeffs.b2;
            coeffs = struct('j_max',[],'b2',[]);
            warning('Concentration overpotential model is based on fuel cell research where the amount of oxygen supplied to the cathode is limiting. Includes experimental numeric parameters that probably cannot be extended for electrolyzers!')
    end
    
    output = struct('name','Ucon','func',Ucon,'coeffs',coeffs);
    
end