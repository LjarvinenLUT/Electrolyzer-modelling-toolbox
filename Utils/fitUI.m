function [fitParams,gof,PlottingCurves] = fitUI(fitFunc,voltage,current,varargin)
% FITUI Fits a given function for a UI curve with given voltage and
% current data to find unknown coefficient values.
%
%   [fitParam,gof,plottingCurves] = FITUI(fitFunc,voltage,current) fits the
%       given function to data using particle swarm optimisation and
%       weighting beginning and end of the dataset.
%       Inputs:
%           fitFunc -- A func object containing a function handle for the
%                   UI curve and all the necessary coefficients, constants
%                   and measured values for the fit.
%           voltage -- Measured voltage values to be used for fitting. Can
%                       be a nx1 column vector containing the data or a nx2 
%                       matrix containing the data in the first column and 
%                       its standard deviation in the second column.
%           current -- Measured current values to be used for fitting. Can
%                       be a nx1 column vector containing the data or a nx2 
%                       matrix containing the data in the first column and 
%                       its standard deviation in the second column.
%
%   [_] = FITUI(_,'method',m) uses the given method for fitting. Available
%       methods are:
%           - 'NLLSE' -- Non-Linear Least Squares Error regression
%           - 'PS' -- Particleswarm optimisation for the least squares
%                       error. Error estimation is performed with
%                       Marchov Chain Monte Carlo -method (MCMC)
%
%   [_] = FITUI(_,'weights',w) allows user to choose whether the beginning
%       and end of the data set are weighted more than the middle. The 
%       reason for the weights is to improve fitting of activation and
%       concentration overpotentials, whose effects are limited to low
%       and high current densities, respectively. Available options
%       are:
%           - 'default' -- Use exponentially changing weight for the 
%                           beginning and end.
%           - 'none' -- Do not use weights for the beginning and end
% 
%   Output:
%       fitParam -- Fit coefficient values in a table with their standard
%                   deviation (confidence value with 1 sigma)
%       gof -- Goodness of fit values in a structure.
%       PlottingCurves -- A structure containing the data used for fitting
%                           as well as curves calculated from the fit.
%       
%   Function FITUI updates the workspace of input func to include the
%   coefficients
%
%	See also FUNC, ELECTROLYZERMODEL


defaultMethod = 'PS';
defaultWeights = 'default';

parser = inputParser;
addRequired(parser,'func',@(x) isa(x,'func'))
addRequired(parser,'voltage',@(x) isnumeric(x))
addRequired(parser,'current',@(x) isnumeric(x))
addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x))
addParameter(parser,'weights',defaultWeights,@(x) ischar(x)||isstring(x))

parse(parser,fitFunc,voltage,current,varargin{:});

method = upper(string(parser.Results.method));
weightsMethod = lower(string(parser.Results.weights));

if length(voltage)~=length(current)
    error('Given voltage and current vectors have incompatible sizes.')
end

if any(current < 0)
    warning("Negative current values may result in poor fit quality if caused by random variation in the data.")
elseif any(current < 0.005)
    warning("Very low current values may result in poor fit quality if caused by random variation in the data.")
end

%% Destructurize function handle, get coefficients and their limits, problem variable names and their values
[funcHandle,coefficients,problemVariableNames,problemVariables] = fitFunc.destructurize('current');

[lb, ub, start] = getArgumentLimits(coefficients, current(:,1));

%% Parse error vectors, if included
if length(voltage(1,:)) == 1 % No measurement standard deviation given
    voltageStd = zeros(size(voltage));
elseif length(voltage(1,:)) == 2 % Measurement standard deviation given as the second column of input matrix
    voltageStd = voltage(:,2);
    voltage = voltage(:,1);
else
    error("Voltage input must be either a column vector or a matrix with two columns containing measured value and its standard deviation")
end

if length(current(1,:)) == 1 % No measurement standard deviation given
    currentStd = zeros(size(current));
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
        weights = (0.01.^x + 0.01.^(1-x) - 0.02)./(1-0.01);
    case "none"
        % Don't apply wieghts
        weights = ones(size(current));
end

% Additinal weights due to measurement accuracy
weights = weights.*errorWeights;

%% Perform fit according to the chosen method	
switch method

    case "NLLSE" % Non-Linear Least Squares Error regression approach

        methodStr = "Non-Linear Least Squares Error Regression";        

        % Retry fitting if R^2 not good enough
        ssr_best = Inf;
        tolX = 1e-8;
        maxIter = 5000;
        for i = 1:10
            fo = fitoptions('Method','NonlinearLeastSquares',...
                'Robust','Bisquare',...
                'Lower',lb,...
                'Upper',ub,...
                'StartPoint',start,...
                'Weights',weights,...
                'MaxIter',maxIter,...
                'TolX', tolX,...
                'Display','off');
            
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
            gofIter = goodnessOfFit(fitVoltage,voltage);
            
            if gofIter.ssr < ssr_best
                ssr_best = gofIter.ssr;
                % Coefficient values
                coeffValues = coeffvalues(fittedCurve);
                gof = gofIter;
                if gofIter.rsquare >0.998 % If good enough fit is achieved:
                    break; % terminate loop
                end
            end
            switch output.exitflag
                case 0 % Fit algorithm exited due to iteration limitation
                    % Induce randomness into starting variables to try find
                    % the global minimum
                    ialpha = strcmp(coefficients,'alpha');
                    start(ialpha) = min(max(rand*(ub(ialpha)-lb(ialpha)),lb(ialpha)),ub(ialpha));
                    ij0 = strcmp(coefficients,'j0');
                    js = logspace(-10,0,1000);
                    start(ij0) = min(max(js(round(rand*1000)),lb(ij0)),ub(ij0));
                case 2 % Fit algorithm exited due to coefficient change limitation
                    tolX = tolX*1e-2; % Increase tolerance
            end
            
        end
        
        % Fit 2 sigma confidence bounds [lower upper]
        sigma = mean(abs(confint(fittedCurve)-coeffValues))/2; 
    
        
    
    case "PS" % Particle swarm approach
        
        methodStr = "Particle Swarm Optimisation";
		
        nvars = length(coefficients); % Amount of fit parameters
        
        % Modify function handle to use vector input for fit parameters
        modFuncHandle = vectorifyFuncHandle(funcHandle,length(coefficients),length(problemVariableNames));
        
        % Creating objective function for particle swarm: sum of square
        % residuals
        fitfun = @(x) sum((modFuncHandle(x,problemVariables{:},current)-voltage).^2.*weights);
        
        % Particle swarm options
        options = optimoptions('particleswarm','SwarmSize',600,'HybridFcn',@fmincon,'Display','off');%,'HybridFcn',@fmincon
        
        % Retry fitting if R^2 not good enough
        fval_best = Inf;
        for i = 1:10
            [coeff,fval,exitflag,output] = particleswarm(fitfun,nvars,lb,ub,options);
            
            % Voltage values obtained from the fit
            fitVoltage = modFuncHandle(coeff,problemVariables{:},current);
            
            % Unweighted goodness of fit values
            gofIter = goodnessOfFit(fitVoltage,voltage);
            
            if fval < fval_best % If previous best fit is improved
                fval_best = fval;
                coeffValues = coeff;
                gof = gofIter;
                if gofIter.rsquare >0.998 % If good enough fit is achieved:
                    break; % terminate loop
                end
            end
        end
        
        %% Call MCMC for uncertainty estimation
        
        coeffStruct = cell(1,length(coefficients));
        for i = 1:length(coefficients)
            coeffStruct{i} = {coefficients{i},coeffValues(i),lb(i),ub(i)};
        end
        
        % Fit 95% confidence bounds [lower upper]
        mcmcfun = @(x) sum((modFuncHandle(x,problemVariables{:},current)-voltage).^2);
        sigma = mcmc(current,voltage,coeffStruct,mcmcfun,gof.ssr);
end

%% Print message about the fit
fprintf('\nData fit performed using %s approach\nR^2: %f\n', methodStr,gof.rsquare)


%% Use map to store fitting parameters can be referenced by name
% Example: fit_param.j0)
fitParams = array2table([coeffValues;sigma], 'VariableNames', coefficients);

%% Input fitted coefficients to the func object
for i = 1:length(coefficients)
    fitFunc.replaceParams(coefficients{i},[coeffValues(i) sigma(i)]);
end

%% Create structure for plotting vectors
fullFitCurrent = min(current):0.001:max(current);
fullFitVoltage = fitFunc.calculate('current',fullFitCurrent);

% Take samples from the dense data vectors
        N = 100; % Number of evenly taken voltage samples
        voltageSamples = linspace(min(fullFitVoltage),max(fullFitVoltage),N)';
%         jsamples = linspace(min(jmeas),jL-0.01,N)'; % Excluding mass transport limitations
%         jmeassamp = nan(N,1);
%         Umeassamp = nan(N,1);
        iii = 1;
        for ii = 1:N
            Udif = abs(fullFitVoltage-voltageSamples(ii));
            [~,ind] = min(Udif);
            if iii == 1 || (iii > 1 && abs(fullFitCurrent(ind)-currentSampled(iii-1))>0.02)
                currentSampled(iii,1) = fullFitCurrent(ind); % Final sampled current vector
                voltageSampled(iii,1) = fullFitVoltage(ind); % Final sampled voltage vector
                iii = iii+1;
            end
        end


PlottingCurves = struct('currentMeasured',[current currentStd],'voltageMeasured',[voltage voltageStd],'currentFit',currentSampled,'voltageFit',voltageSampled);

end


%%
function [lower, upper, start] = getArgumentLimits(argumentList, I)
% GETARGUMENTLIMITS gets lower and upper limits and also startpoint of each 
%	fitting coefficient.
%       Inputs:
%           argumentList -- List of coefficient names in a cell array.
%           I -- Current density for determining the limits for limiting
%                   current density, 'j_lim'.
%       Outputs:
%           lower -- Array of the lower limits for the coefficients in the
%                       same order as listed in argumentList.
%           higher -- Array of the higher limits for the coefficients in 
%                       the same order as listed in argumentList.
%           start -- Array of starting points for non-linear least squares
%                       error estimation listed in the same order as in
%                       argumentsList

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
function modFuncHandle = vectorifyFuncHandle(funcHandle,nCoeffs,nProbVars)
% VECTORIFYFUNCHANDLE returns a function handle with one vector for the 
%   coefficients transformed from input function handle with separate 
%   coefficients that are organized in order 
%   (coefficients, problem variables, independent variable).
%       Inputs:
%           funcHandle -- The function handle to modify
%           nCoeffs -- Number of coefficients for the fit function
%           nProbVars -- Number of problem variables
%       Output:
%           modFuncHandle -- Function handle modified to take cell array
%                               inputs for coefficients instead of listing 
%                               them separately.
       
    % Creating symbolic variables for the function handle
    x = num2cell(sym('x', [1 nCoeffs])); % fit parameters
    y = num2cell(sym('y', [1 nProbVars])); % problem variables
    syms current; % Current input
    
    % Calculating the symbolic result
    z = funcHandle(x{:},y{:},current);
    
    % Creating a modified function handle from the symbolic result
    modFuncHandle = matlabFunction(z,'Vars',[{cell2sym(x)},y(:)',{current}]); 
    
end


%%
function gof = goodnessOfFit(yfit,y)
% GOODNESSOFFIT returns a structure containing unweighted:
%   - SSR (Sum of Square Error)
%   - RMSE (Root Mean Square Error)
%   - R^2 value
%   of the fit
%       Inputs:
%           yfit -- y values calculated based on the fit.
%           y -- Measured y values.
%       Output:
%           gof -- A structure containing the abovementioned goodness of
%           fit values.

residuals = y-yfit; % Residuals
divFromMean = y-mean(y); % Divergence from mean

ssr = sum(residuals.^2); % Sum of square residuals
rmse = sqrt(mean(residuals.^2)); % Root mean square error
rsquare = 1 - ssr/sum(divFromMean.^2); % R^2 value

gof = struct('ssr',ssr,'rmse',rmse,'rsquare',rsquare); % Goodness of fit structure
end


%%
function sigma = mcmc(x,y,coeffs,ssfun,ssr)
% MCMC estimates the posterior distribution for fit coefficients using
%   Markov Chain Monte Carlo
%
% the coeffs structure: name, init.val, min, max
% coeffs = {
% {'\theta_1',initVal,lowB,upB}
% {'\theta_2',initVal,lowB,upB}
% };

n = length(x); % Length of the data
data = [x,y]; % Data matrix

% Model
% if length(coeffs)
%     
% end
Model.ssfun = ssfun; % Sum of squares function
Model.sigma2 = ssr/(n-2); % Proposal variance

% Options
Options.nsimu = 20000; % number of samples
Options.qcov = 0.01*eye(n); % (initial) proposal covariance
Options.method = 'dram'; % method: DRAM
Options.adaptint = 100; % adaptation interval
Options.verbosity = 0; % Suppress output

% calling MCMCRUN
[results,chain] = mcmcrun(Model,data,coeffs,Options);

% visualizing the results using MCMCPLOT
% figure
% mcmcplot(chain,[],results.names);
% figure
% mcmcplot(chain,[],results.names,'pairs');

y = chainstats(chain,results,0);

sigma = y(:,2)';

end
