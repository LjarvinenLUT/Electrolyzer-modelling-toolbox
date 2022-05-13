function [fitParams,gof] = fitUI(fitFunc,voltage,current,varargin)
% FITUI Fits a given function for a UI curve with given voltage and
% current data to find unknown parameter values.
%
%   [fitParam,gof] = FITUI(fitFunc,voltage,current) fits the
%       given function to data using particle swarm optimisation and
%       weighting beginning and end of the dataset.
%       Inputs:
%           fitFunc -- A func object containing a function handle for the
%                   UI curve and all the necessary parameters, constants
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
%       and high current densities, respectively. Equation used for
%       weighting is:
%           weight = tan(y).^2 + 0.1,
%       where y is determined differently based on the user choise.
%       
%       Available options are:
%           - 'hl' or 'lh' -- For weighting both low and high current
%                               values. y is the current normalized to
%                               interval [-pi/4 pi/4]
%           - 'h' -- For weighting high current values. y is the curent
%                       normalized to interval [0 pi/4]
%           - 'l' -- For weighting low current values. y is the curent
%                       normalized to interval [-pi/4 0]
%           - 'none' -- No additional weights applied (default)
%
%   Output:
%       fitParam -- Fit parameter values in a table with their standard
%                   deviation (confidence value with 1 sigma)
%       gof -- Goodness of fit values in a structure.
%
%   Function FITUI updates the workspace of input func to include the
%   parameters
%
%	See also FUNC, ELECTROLYZERMODEL


defaultMethod = 'PS';
defaultWeights = 'none';

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

if ~strcmpi(method, 'PS') && ~strcmpi(method, 'NLLSE')
    error("Given fitting method: " + method + " is not available. The available methods are PS (particle swarm) and NLLSE (non-linear least squares error).")
end



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
    errorWeights = min((1./std.^2)/1000,1);
end


%% Weighting
addedWeights = addWeight(current,weightsMethod);

% Additinal weights due to measurement accuracy
weights = addedWeights.*errorWeights;


%% Perform fit according to the chosen method
switch method
    
    case "NLLSE" % Non-Linear Least Squares Error regression approach
        
        methodStr = "Non-Linear Least Squares Error Regression";
        
        %Destructurize function handle, get parameters and their limits, problem variable names and their values
        [funcHandle,parameters,problemVariableNames,problemVariables] = fitFunc.destructurize('current');
        
        [lb, ub, start] = getArgumentLimits(parameters,...
            fitFunc.Fitlims,...
            current);
        
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
            
            fitFun = fittype(funcHandle,...
                'dependent','voltage',...
                'coefficients',parameters,...
                'independent','current',...
                'problem',problemVariableNames,...
                'options',fo);
            
            [fittedCurve,~,output] = fit(current,...
                voltage,...
                fitFun,...
                'problem',problemVariables);
            
            
            % Voltage values obtained from the fit
            fitVoltage = fittedCurve(current);
            
            % Unweighted goodness of fit values
            gofIter = goodnessOfFit(fitVoltage,voltage);
            
            if gofIter.ssr < ssr_best
                ssr_best = gofIter.ssr;
                % Parameter values
                paramValues = coeffvalues(fittedCurve);
                gof = gofIter;
                if gofIter.rsquare >0.998 % If good enough fit is achieved:
                    break; % terminate loop
                end
            end
            switch output.exitflag
                case 0 % Fit algorithm exited due to iteration limitation
                    % Induce randomness into starting variables to try find
                    % the global minimum
                    ialpha = strcmp(parameters,'alpha');
                    start(ialpha) = min(max(rand*(ub(ialpha)-lb(ialpha)),lb(ialpha)),ub(ialpha));
                    ij0 = strcmp(parameters,'j0');
                    js = logspace(-10,0,1000);
                    start(ij0) = min(max(js(round(rand*1000)),lb(ij0)),ub(ij0));
                case 2 % Fit algorithm exited due to parameter change limitation
                    tolX = tolX*1e-2; % Increase tolerance
            end
            
        end
        
        % 2-sigma confidence bounds
        sigma = mean(abs(confint(fittedCurve)-paramValues));
        
        
        
    case "PS" % Particle swarm approach
        
        methodStr = "Particle Swarm Optimisation";
        
        %% Vectorify function handle, get parameters and their limits, problem variable names and their values
        [funcHandle,parameters,~,problemVariables] = fitFunc.vectorify('current');
        
        [lb, ub, ~] = getArgumentLimits(parameters,...
            fitFunc.Fitlims,...
            current);
        
        nvars = length(parameters); % Amount of fit parameters
        
        % Creating objective function for particle swarm: sum of square
        % residuals
        fitFun = @(x) sum((funcHandle(x,problemVariables{:},current)-voltage).^2.*weights);
        
        % Particle swarm options
        options = optimoptions('particleswarm',...
            'SwarmSize',600,...
            'HybridFcn',@fmincon,...
            'Display','off');%,'HybridFcn',@fmincon
        
        % Retry fitting if R^2 not good enough
        fval_best = Inf;
        for i = 1:10
            [param,fval,exitflag,output] = particleswarm(fitFun,...
                nvars,...
                lb,...
                ub,...
                options);
            
            % Voltage values obtained from the fit
            fitVoltage = funcHandle(param,problemVariables{:},current);
            
            % Unweighted goodness of fit values
            gofIter = goodnessOfFit(fitVoltage,voltage);
            
            if fval < fval_best % If previous best fit is improved
                fval_best = fval;
                paramValues = param;
                gof = gofIter;
                if gofIter.rsquare >0.998 % If good enough fit is achieved:
                    break; % terminate loop
                end
            end
        end
        
        %% Call MCMC for uncertainty estimation
        
        params = cell(1,length(parameters));
        for i = 1:length(parameters)
            params{i} = {parameters{i},paramValues(i),lb(i),ub(i)};
        end

        % If the amount of measurements is too low and measurement error is
        % significant, increase measurement points artificially from the
        % measurement error (gaussian distribution around existing
        % measurements) DISABLED
        if false && length(current) < 20 && ...
                (any(abs(currentStd-current)>1e-5) ||...
                any(abs(voltageStd-voltage)>1e-5))
            [currentNew,voltageNew] = increaseNumberOfPoints(current,...
                currentStd,...
                voltage,...
                voltageStd,...
                2*length(current));
        else
            currentNew = current;
            voltageNew = voltage;
        end

        % Fit 95% confidence bounds
        mcmcfun = @(x) sum((funcHandle(x,problemVariables{:},currentNew)-voltageNew).^2);
        [results,chain] = mcmc(currentNew,voltageNew,params,mcmcfun,gof.ssr);

%         % visualizing the results using MCMCPLOT
%         % TODO add user choice
%         figure('Name','Covariance')
%         mcmcplot(chain,[],results,'pairs');
% 
%         figure('Name','Chains')
%         subplot(2,2,1)
%         plot(chain(:,1))
%         legend(params{1}{1})
%         subplot(2,2,2)
%         plot(chain(:,2))
%         legend(params{2}{1})
%         subplot(2,2,3)
%         plot(chain(:,3))
%         legend(params{3}{1})
%         subplot(2,2,4)
%         plot(chain(:,4))
%         legend(params{4}{1})
% 
%         figure('Name','Distribution')
%         mcmcplot(chain,[],results,'denspanel',4);

%         model = @(x) (funcHandle(x,problemVariables{:},currentNew));
%         out = mcmcpred(results,chain(50:100:end,:),[],current,model);
%         figure('Name','Plot with distribution')
%         mcmcpredplot(out);

        y = chainstats(chain,results,0); % [value, std, something, something]
        
        % 2-sigma confidence bounds
        sigma = y(:,2)'*2;
end

%% Print message about the fit
fprintf('\nData fit performed using %s approach\nR^2: %f\n', methodStr,gof.rsquare)


%% Use map to store fitting parameters can be referenced by name
% Example: fit_param.j0)
fitParams = array2table([paramValues;sigma], 'VariableNames', parameters);

%% Input fitted parameters to the func object
for i = 1:length(parameters)
    fitFunc.replaceInWorkspace(parameters{i},[paramValues(i) sigma(i)],'rebuild');
end



end


%%
function [lower, upper, start] = getArgumentLimits(argumentList,LimsStruct, x)
% GETARGUMENTLIMITS gets lower and upper limits and also startpoint of each
%	fitting parameter.
%       Inputs:
%           argumentList -- List of parameter names in a cell array.
%           LimsStruct -- The structure for the limits in the func object
%           x -- Current density for determining the limits for
%                   parameters with limits depending on the independent
%                   variable.
%       Outputs:
%           lower -- Array of the lower limits for the parameters in the
%                       same order as listed in argumentList.
%           higher -- Array of the higher limits for the parameters in
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
    try
        limsList = LimsStruct.(argumentList{i});
        lims = nan(1,3);
        for j = 1:3
            if isnumeric(limsList{j})
                lims(j) = limsList{j};
            elseif ischar(limsList{j})||isstring(limsList{j})
                try
                    lims(j) = eval(limsList{j});
                catch ME
                    if strcmp(ME.identifier,'MATLAB:UndefinedFunction')
                        errormsg = "Expression for the fit limit of "+...
                            argumentList{i} + " contains unrecognized "+...
                            "functions or variables. Please recheck "+...
                            "the expression and ensure that the dependent "+...
                            "variable is marked with x.";
                        error(errormsg);
                    else
                        rethrow(ME)
                    end
                end
            else
                error("Fit limits for " + argumentList{i} + " are of unsupported type.")
            end
        end
                
        lower(i) = lims(1);
        upper(i) = lims(3);
        start(i) = lims(2);
        
    catch ME
        if strcmp(ME.identifier,'MATLAB:nonExistentField')
            errormsg = "fit_UI.getArgumentLimits: No parameter limits "...
                +"could be found for "+ argumentList{i} +". Consider "...
                +"setting limits for all the parameters using method "...
                +"func.setFitlims.";
            error(errormsg);
        else
            rethrow(ME)
        end
    end
        
end
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
function [results,chain] = mcmc(x,y,params,ssfun,ssr)
% MCMC estimates the posterior distribution for fit parameters using
%   Markov Chain Monte Carlo
%
% the params structure: name, init.val, min, max
% params = {
% {'\theta_1',initVal,lowB,upB}
% {'\theta_2',initVal,lowB,upB}
% };

nSamples = length(x); % Length of the data
nParams = length(params); % Number of parameters
dof = nSamples-nParams; % Degrees of freedome
data = [x,y]; % Data matrix

% Setup Model
Model.ssfun = ssfun; % Sum of squares function
Model.sigma2 = ssr/dof; % Proposal variance

    
% Options
Options.nsimu = nSamples*1000; % number of samples
Options.qcov = 0.01*eye(nSamples); % (initial) proposal covariance
Options.method = 'dram'; % method: DRAM
Options.adaptint = 100; % adaptation interval
Options.verbosity = 0; % Suppress output
Options.waitbar = 0; % Suppress waitbar

% Run Metropolis-Hastings MCMC simulation
[results,chain] = mcmcrun(Model,data,params,Options);

end


%%
function [xNew,yNew] = increaseNumberOfPoints(x,xStd,y,yStd,n)
% INCREASENUMBEROFPOINTS artificially increases the number of measurement
% points based on measurement error. Assumes gaussian distribution around
% the measured value. Generates points randomly around the real
% measurements. 

xNew = x;
yNew = y;

if any(length(x)~=[length(y) length(xStd) length(yStd)])
    error("Input vector lengths are not equal")
elseif length(x) >= n
    return;
end

ki = length(x);
k = ki;

while k < n
    index = randi(ki);
    xNew(k+1) = mean(x(index)+randn(25,1)*xStd(index)); 
    yNew(k+1) = mean(y(index)+randn(25,1)*yStd(index));
    k = k+1;
end

end

%%
function addedWeights = addWeight(current,weightsMethod)
% ADDWEIGHT calculates the weighting curve used in the fitting
switch weightsMethod
    case {"hl","lh"}
        % Weigh beginning and end of the measured current spectrum.
        % y = current normalized to the interval [-pi/4 pi/4]
        x = (current-(max(current)+min(current))/2).^2;
    case "h"
        % Weigh the beginning of the measured current spectrum
        % y = current normalized to the interval [0 pi/4]
        x = (current-min(current)).^4;
    case "l"
        % Weigh the end of the measured current spectrum
        % y = current normalized to the interval [-pi/4 0]
        x = (current-max(current)).^4;
    case "none"
        % Don't apply weights
        x = ones(size(current));
    otherwise
        error("Nonvalid weight method! Set the method as one of the following: 'l' (low), 'h' (high), 'hl'/'lh' (high and low) or 'none'")
end

y = x/max(abs(x))*pi/4;
addedWeights = tan(y).^2 + 0.1;

end