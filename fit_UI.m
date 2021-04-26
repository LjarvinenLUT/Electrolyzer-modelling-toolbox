
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


%% Perform fit according to the chosen method	
switch method
    
    case "NLLSE" % Non-Linear Least Squares Error regression approach
        % weight beginning and end
        x = (0:length(Current)-1)';
        weights = exp(x).^(-7/length(x)) + exp(x-x(end)).^(7/length(x));

        
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
        
        [fitted_curve,gof,output] = fit(Current,Voltage,fitfun);
        
        coeff_values = coeffvalues(fitted_curve);
        err_bounds = confint(fitted_curve)-coeff_values; % Fit 95% confidence bounds [lower upper]
      
    case "PS" % Particle swarm approach
		
        nvars = length(coefficients); % Amount of fit parameters
        
        % Modify function handle to use vector input for fit parameters
        mod_func_handle = vectorify_func_handle(func_handle,nvars);
        
        % Creating objective function for particle swarm: sum of square
        % residuals
        fitfun = @(x) sum((mod_func_handle(x,Current)-Voltage).^2);
        
        % Particle swarm options
        options = optimoptions('particleswarm','SwarmSize',600,'HybridFcn',@fmincon);%,'HybridFcn',@fmincon
        
        % Minimisation of the objective function using particle swarm: 
        % best of three
        fval_best = Inf;
        for i = 1:3
            [coeff,fval,exitflag,output] = particleswarm(fitfun,nvars,lb,ub,options);
            if fval < fval_best
                fval_best = fval;
                coeff_values = coeff;
            end
        end
        
        % Goodness of fit values
        sse = fval_best;
        rmse = sqrt(mean((mod_func_handle(coeff_values,Current)-Voltage).^2));
        rsquare = 1 - sum(((Voltage-mod_func_handle(coeff_values,Current)).^2).^2)/sum((Voltage-mean(Voltage)).^2);
        
        err_bounds = nan(size(coeff_values));
        gof = [];
end


% Use map to store fitting parameters can be referenced by name
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
