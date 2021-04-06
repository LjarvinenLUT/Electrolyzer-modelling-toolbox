
function fit_param = fit_UI(func_handle,Voltage,Current,varargin)

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
[lower, upper, start] = getArgumentLimits(coefficients, Voltage, Current);
	
	
switch method
    
    case "NLLSE" % Non-Linear Least Squares Error regression approach
        % weight beginning and end
        x = (0:length(Current)-1)';
        weights = exp(x).^(-5/length(x)) + exp(x-x(end)).^(5/length(x));

        
        fo = fitoptions('Method','NonlinearLeastSquares',...
			'Lower',lower,...
			'Upper',upper,...
			'StartPoint',start,...
            'Weights',weights,...
            'MaxIter',1500,...
            'Display','notify');
        
        fitfun = fittype(func_handle,...
			'dependent','Voltage',...
			'coefficients',coefficients,...
			'independent','Current',...
			'options',fo);
        
        [fitted_curve,gof] = fit(Current,Voltage,fitfun);
        
        fit_param = coeffvalues(fitted_curve);
        
      
    case "PS" % Particle swarm approach

        % Not functional
		
        fun = @(x) sum((func_handle(x(1),x(2),x(3),x(4),Current)-Voltage).^2);
        nvars = 4;
        lb = [1e-10,0,0,max(Current)];
        ub = [1,1,40,5];
        
        options = optimoptions('particleswarm','SwarmSize',200,'HybridFcn',@fmincon);
        
        [fit_param,fval,exitflag,output] = particleswarm(fun,nvars,lb,ub,options);
        
end


% Use map to store fitting parameters can be referenced by name
% Example: fit_param('j0')
fit_param = containers.Map(coefficients, coeff_values);
end


