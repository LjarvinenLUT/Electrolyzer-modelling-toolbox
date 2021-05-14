
function [fit_param,err_bounds,gof] = fit_UI(func,voltage,current,varargin)

addpath Utils;
% Add path to toolbox with MCMC
addpath Utils\mcmcstat;

defaultMethod = 'PS';
defaultWeights = 'default';

parser = inputParser;
addRequired(parser,'func',@(x) isa(x,'func'))
addRequired(parser,'voltage',@(x) isnumeric(x))
addRequired(parser,'current',@(x) isnumeric(x))
addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x))
addParameter(parser,'weights',defaultWeights,@(x) ischar(x)||isstring(x))

parse(parser,func,voltage,current,varargin{:});

method = upper(string(parser.Results.method));
weightsMethod = lower(string(parser.Results.weights));


%% Destructurize function handle, get coefficients and their limits, problem variable names and their values
[funcHandle,coefficients,problemVariableNames,problemVariables] = func.destructurize('current');

[lb, ub, start] = getArgumentLimits(coefficients, current);

%% Parse error vectors, if included
if length(voltage(1,:)) == 1 % No measurement standard deviation given
    voltageStd = 0;
elseif length(voltage(1,:)) == 2 % Measurement standard deviation given as the second column of input matrix
    voltageStd = voltage(:,2);
    voltage = voltage(:,1);
else
    error("Voltage input must be either a column vector or a matrix with two columns containing measured value and its standard deviation")
end

if length(current(1,:)) == 1 % No measurement standard deviation given
    currentStd = 0;
elseif length(current(1,:)) == 2 % Measurement standard deviation given as the second column of input matrix
    currentStd = current(:,2);
    current = current(:,1);
else
    error("Current input must be either a column vector or a matrix with two columns containing measured value and its standard deviation")
end

% Total error obtained from the squared sum
std = sqrt(voltageStd.^2 + currentStd.^2);
if std == 0
    errorWeights = 1;
else
    errorWeights = 1./std.^2;
end


%% Weighting
switch weightsMethod
    case "default"
        % Weigh beginning and end of the measured current spectrum
        x = current/max(current);
        weights = (exp(log(0.1^2)*x) + exp(-log(0.1^2)*(x-x(end))) - 2*exp(log(0.1^2)))./(1-exp(log(0.1^2)));
    case "none"
        % Don't apply wieghts
        weights = ones(size(current));
end

% Additinal weights due to measurement accuracy
weights = weights.*errorWeights;

%% Perform fit according to the chosen method	
switch method

    case "NLLSE" % Non-Linear Least Squares Error regression approach

        method_str = "Non-Linear Least Squares Error Regression";
        
        fo = fitoptions('Method','NonlinearLeastSquares',...
			'Lower',lb,...
			'Upper',ub,...
			'StartPoint',start,...
            'Weights',weights,...
            'MaxIter',1500,...
            'Display','notify');
        
        fitfun = fittype(funcHandle,...
			'dependent','voltage',...
			'coefficients',coefficients,...
			'independent','current',...
            'problem',problemVariableNames,...
			'options',fo);
        
        [fittedCurve,~,output] = fit(current,voltage,fitfun,'problem',problemVariables);

        
        % Voltage values obtained from the fit
        fitVoltage = fittedCurve(current);
        
        % Unweighted goodness of fit values
        gof = goodness_of_fit(fitVoltage,voltage);
        
        % Coefficient values
        coeffValues = coeffvalues(fittedCurve);       
        
        % Fit 95% confidence bounds [lower upper]
        err_bounds = confint(fittedCurve)-coeffValues; 
    
        
    
    case "PS" % Particle swarm approach
        
        method_str = "Particle Swarm Optimisation";
		
        nvars = length(coefficients); % Amount of fit parameters
        
        % Modify function handle to use vector input for fit parameters
        modFuncHandle = vectorifyFuncHandle(funcHandle,coefficients,problemVariableNames);
        
        % Creating objective function for particle swarm: sum of square
        % residuals
        fitfun = @(x) sum((modFuncHandle(x,problemVariables,current)-voltage).^2.*weights);
        
        % Particle swarm options
        options = optimoptions('particleswarm','SwarmSize',600,'HybridFcn',@fmincon,'Display','off');%,'HybridFcn',@fmincon
        
        % Retry fitting if R^2 not good enough
        fval_best = Inf;
        for i = 1:5
            [coeff,fval,exitflag,output] = particleswarm(fitfun,nvars,lb,ub,options);
            
            % Voltage values obtained from the fit
            fitVoltage = modFuncHandle(coeff,problemVariables,current);
            
            % Unweighted goodness of fit values
            gof_iter = goodness_of_fit(fitVoltage,voltage);
            
            if gof_iter.rsquare >0.995 % If good enough fit is achieved:
                coeffValues = coeff;
                gof = gof_iter;
                break; % terminate loop
            elseif fval < fval_best % If previous best fit is improved
                fval_best = fval;
                coeffValues = coeff;
                gof = gof_iter;
            end
            
        end
        
        % Fit 95% confidence bounds [lower upper] (TODO)
        err_bounds = nan(size(coeffValues));
end

%% Print message about the fit
fprintf('\nData fit performed using %s approach\nR^2: %f\n', method_str,gof.rsquare)


%% Use map to store fitting parameters can be referenced by name
% Example: fit_param.j0)
fit_param = array2table(coeffValues, 'VariableNames', coefficients);

% Input fitted coefficients to the func object
for i = 1:length(coefficients)
    func.Workspace.Coefficients.(coefficients{i}) = coeffValues(i);
end
end


%%
function [lower, upper, start] = getArgumentLimits(argumentList, I)
%GETARGUMENTLIMITS gets lower and upper limits and also startpoint of each 
%fitting parameter

% Initialize arrays for lower and upper limits and for start points
lower = zeros(1, length(argumentList));
upper = zeros(1, length(argumentList));
start = zeros(1, length(argumentList));

% Get limits for each function argument if none is found give an error
for i = 1:length(argumentList)
    switch argumentList{i}
        case 'j0'
            lower(i) = 1e-10;
            upper(i) = 1;
            start(i) = 1e-5;
        case 'alpha'
            lower(i) = 0;
            upper(i) = 1;
            start(i) = 0.1;
        case 'r'
            lower(i) = 0;
            upper(i) = inf;
            start(i) = 1;
        case 'j_lim'
            lower(i) = max(I);
            upper(i) = 3*max(I);
            start(i) = max(I)*1.01;
        case 'Uerr'
            lower(i) = -inf;
            upper(i) = inf;
            start(i) = 0;
        otherwise
            error(['fit_UI.getArgumentLimits: No argument limits could be found for ' argumentList{i}]);
    end
end
end




%% 
function mod_func_handle = vectorifyFuncHandle(func_handle,coeffs,probVarsNames)
% VECTORIFY_FUNC_HANDLE returns a function handle with one vector for the 
% coefficients transformed from input function handle with separate 
% coefficients 
       
    % Creating symbolic variables for the function handle
    x = num2cell(sym('x', [1 length(coeffs)])); % fit parameters
    y = num2cell(sym('y', [1 length(probVarsNames)])); % problem variables
    syms current; % Current input
    
    % Calculating the symbolic result
    z = func_handle(x{:},y{:},current);
    
    % Creating a modified function handle from the symbolic result
    mod_func_handle = matlabFunction(z,'Vars',{cell2sym(x),cell2sym(y),current});
    
end

%%
function gof = goodness_of_fit(yfit,y)
% GOODNESS_OF_FIT returns a structure containing unweighted:
% - SSR (Sum of Square Error)
% - RMSE (Root Mean Square Error)
% - R^2 value
% of the fit

residuals = y-yfit; % Residuals

sse = sum(residuals.^2); % Sum of square residuals
rmse = sqrt(mean(residuals.^2)); % Root mean square error
rsquare = 1 - sum(residuals.^2)/sum((y-mean(y)).^2); % R^2 value

gof = struct('sse',sse,'rmse',rmse,'rsquare',rsquare); % Goodness of fit structure
end
