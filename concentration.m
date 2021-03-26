% Concentration overpotetial function handle generator
% Inputs:       model - Used model reference, numeric, from which article
%               p_O2 - Oxygen partial pressure
%               p_sat - Vapor saturation pressure
%               T - Temperature

function Ucon = concentration(varargin)
    
    defaultModel = 1;

    parser = inputParser;
    addRequired(parser,'p_O2',@(x) isnumeric(x))
    addRequired(parser,'p_sat',@(x) isnumeric(x))
    addRequired(parser,'T',@(x) isnumeric(x))
    addParameter(parser,'model',defaultModel,@(x) isnumeric(x)&&isscalar(x))
    
    parse(parser,varargin{:});
    
    p_O2 = parser.Results.p_O2;
    p_sat = parser.Results.p_sat;
    T = parser.Results.T;
    model = parser.Results.model;
    
    %% Global parameters
    [F,R,n_e] = get_constants();

    f = R/(n_e*F);
    
    %% Error checking
    
    
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
            Ucon = @(Cx1,Cx0) f*T*log(Cx1/Cx0);
        case 2 % ???
            Ucon = @(j,b,j_lim) (R*T)/(2*b*F)*log(1 + j/j_lim);
        case 3 % J. Pukrushpan "Modeling and control for PEM fuel cell stack system" https://doi.org/10.1109/ACC.2002.1025268
            px = p_O2/0.1173 + p_sat;
            atm = 1.01325;
            if px < 2 * atm
                b1 = (7.16*10^-4*T - 0.622)*px + (-1.45*10^-3*T + 1.68); 
            else
                b1 = (8.66*10^-5*T - 0.068)*px + (-1.6*10^-4*T + 0.54);
            end
            
            Ucon = @(j,b2,j_max) j*(b1*j/j_max)^b2;
    end
end