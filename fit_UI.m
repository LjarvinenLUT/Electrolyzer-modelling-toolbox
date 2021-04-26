
function [fit_param,err_bounds,gof] = fit_UI(func_handle,Voltage,Current,varargin)

defaultMethod = 'PS';

parser = inputParser;
addRequired(parser,'func_handle',@(x) isa(x,'function_handle'))
addRequired(parser,'Voltage',@(x) isnumeric(x))
addRequired(parser,'Current',@(x) isnumeric(x))
addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x))

parse(parser,func_handle,Voltage,Current,varargin{:});

method = upper(string(parser.Results.method));


%% Get coefficients and their limits
coefficients = getFunctionArguments(func_handle);
[lb, ub, start] = getArgumentLimits(coefficients, Voltage, Current);

%% Weigh beginning and end of the measured current spectrum
x = Current/max(Current);
weights = exp(log(0.1^2)*x) + exp(-log(0.1^2)*(x-x(end)));

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
        
        fitfun = fittype(func_handle,...
			'dependent','Voltage',...
			'coefficients',coefficients,...
			'independent','Current',...
			'options',fo);
        
        [fitted_curve,~,output] = fit(Current,Voltage,fitfun);

        
        % Voltage values obtained from the fit
        fit_Voltage = fitted_curve(Current);
        
        % Unweighted goodness of fit values
        gof = goodness_of_fit(fit_Voltage,Voltage);
        
        % Coefficient values
        coeff_values = coeffvalues(fitted_curve);       
        
        % Fit 95% confidence bounds [lower upper]
        err_bounds = confint(fitted_curve)-coeff_values; 
    
        
        
    case "PS" % Particle swarm approach
        
        method_str = "Particle Swarm Optimisation";
		
        nvars = length(coefficients); % Amount of fit parameters
        
        % Modify function handle to use vector input for fit parameters
        mod_func_handle = vectorify_func_handle(func_handle,nvars);
        
        % weight beginning and end
        x = Current/max(Current);
        weights = exp(log(0.1^2)*x) + exp(-log(0.1^2)*(x-x(end)));
        
        % Creating objective function for particle swarm: sum of square
        % residuals
        fitfun = @(x) sum((mod_func_handle(x,Current)-Voltage).^2.*weights);
        
        % Particle swarm options
        options = optimoptions('particleswarm','SwarmSize',600,'HybridFcn',@fmincon,'Display','off');%,'HybridFcn',@fmincon
        
        % Retry fitting if R^2 not good enough
        fval_best = Inf;
        for i = 1:5
            [coeff,fval,exitflag,output] = particleswarm(fitfun,nvars,lb,ub,options);
            
            % Voltage values obtained from the fit
            fit_Voltage = mod_func_handle(coeff,Current);
            
            % Unweighted goodness of fit values
            gof_iter = goodness_of_fit(fit_Voltage,Voltage);
            
            if gof_iter.rsquare >0.995 % If good enough fit is achieved:
                coeff_values = coeff;
                gof = gof_iter;
                break; % terminate loop
            elseif fval < fval_best % If previous best fit is improved
                fval_best = fval;
                coeff_values = coeff;
                gof = gof_iter;
            end
            
        end
        
        % Fit 95% confidence bounds [lower upper] (TODO)
        err_bounds = nan(size(coeff_values));
end

%% Print message about the fit
fprintf('Data fit performed using %s approach\nR^2: %f\n', method_str,gof.rsquare)


%% Use map to store fitting parameters can be referenced by name
% Example: fit_param.j0)
fit_param = array2table(coeff_values, 'VariableNames', coefficients);
end

%%
function [coefficients] = getFunctionArguments(theFcn)
    %GETFUNCTIONARGUMENTS returns the fucntion arguments from the provided
    %function as a cell array.
    
    % Get the amount of arguments
    numArguments = nargin( theFcn );

    % Get the string description of the function
    functionString = func2str( theFcn );

    % Allocate space for the cell-string. One of the arguments is the
    % independent fit variable which is not needed here
    coefficients = cell( 1, numArguments - 1);

    % The plan is to move a pair of indices along the "function string" field
    % looking for the commas. The names we want will be between these indices.
    %
    % We know from the form of the string that the first name starts two
    % characters after the "@"
    indexOfAtSign = find( functionString == '@', 1, 'first' );
    ai = indexOfAtSign + 2;
    % Therefore the first comma must be no sooner that the fourth character
    bi = 4;
    % When we start we have found no arguments
    numFound = 0;
    % We will keep looping until we have found all the arguments we expect
    while numFound < numArguments
        % If we have found the end of a argument name
        if functionString(bi) == ',' || functionString(bi) == ')'
            % then increment the "numFound" counter
            numFound = numFound+1;
            % ... store the name
            if ~strcmpi(functionString(ai:(bi-1)), "Current")
                coefficients{numFound} = functionString(ai:(bi-1));
            end
            % and increment the start index
            ai = bi+1;
            % Since the end must be beyond the start, we set the end index
            % beyond the start index.
            bi = ai+1;
        else
            % Otherwise increment the end index
            bi = bi+1;
        end
    end
end

%%
function [lower, upper, start] = getArgumentLimits(argumentList, U, I)
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
        case 'jL'
            lower(i) = max(I);
            upper(i) = 3*max(I);
            start(i) = max(I)*1.01;
        case 'Uerr'
            lower(i) = -inf;
            upper(i) = inf;
            start(i) = 0;
        otherwise
            error('fit_UI.getArgumentLimits: No argument limits could be found for ' + argumentList{i});
    end
end
end




%% 
function mod_func_handle = vectorify_func_handle(func_handle,nvars)
% VECTORIFY_FUNC_HANDLE returns a function handle with one vector for the 
% coefficients transformed from input function handle with separate 
% coefficients 
       
    % Creating symbolic variables for the function handle
    x = num2cell(sym('x', [1 nvars])); % fit parameters
    syms Current; % Current input
    
    % Calculating the symbolic result
    z = func_handle(x{:},Current);
    
    % Creating a modified function handle from the symbolic result
    mod_func_handle = matlabFunction(z,'Vars',{cell2sym(x),Current});
    
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
