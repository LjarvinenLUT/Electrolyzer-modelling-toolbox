function [SynData,FullData,Workspace] = createSyntheticUI(varargin)
% CREATESYNTHETICUI Creates synthetic Voltage-Current pairs from an
%   electrolyzer UI curve. The function simulates a true measurement and
%   thus enables the user to define measurement errors for both current and
%   voltage, as well as the number of independent measurements for each
%   data point to be averaged.
%
%   [SynData,FullData,Workspace] = CREATESYNTHETICUI() creates synthetic UI curve
%       with N = 20 samples using default coefficients. Buttler-Volmer
%       equation is used for activation overpotential. The synthetic data
%       is outputed in SynData structure containing fields 'current' and
%       'voltage', both of them being Nx2 matrices (first column = data,
%       second column = standard deviation). Workspace structure includes
%       the variables, constants and coefficients used for creating the
%       synthetic data.
%
%   [_] = CREATESYNTHETICUI(N) creates synthetic UI curve with user defined
%       number of samples, N.
%
%   [_] = CREATESYNTHETICUI(_,name,value) creates synthetic UI curve using
%       the user defined parameters inputed as name--value pairs. The
%       parameters that are recognised by the function are listed below:
%           - 'eModel', enables user defined model and variables to be used
%               for creating the synthetic data.
%           - 'jLims', enables defining of current limits for the synthetic
%               UI data. Input should be in a form of 1x2 vector 
%               [iLow iHigh].
%               Default: [0.02 5]
%           - 'jErr', enables defining of measurement error for single
%               current measurements. Provided as fraction from the
%               measured value.
%               Default: 0
%           - 'uLims', enables defining of voltage limits for the synthetic
%               UI data. Input should be in a form of 1x2 vector 
%               [uLow uHigh].
%               Default: [0 inf]
%           - 'uErr', enables defining of measurement error for single
%               voltage measurements. Provided as fraction from the
%               measured value.
%               Default: 0
%           - 'nSamp', number of "measurements" from normal distribution
%               per data sample. 
%               Default: 200
%
%   Output:
%       - SynData -- Structure containing sampled data in fields 'current' 
%           and 'voltage', both of them being Nx2 matrices:
%           (first column = data, second column = standard deviation).
%       - FullData -- Structure containing non-sampled data in fields
%           'current' and 'voltage', both of them being Nx2 matrices:
%           (first column = data, second column = standard deviation).
%       - Workspace -- Structure including the variables, constants and
%           coefficients used for creating the synthetic data.
%
%   See also FITUI, FUNC, ELECTROLYZERMODEL

fprintf("------------------------------------------------------------\n"...
    + "Creating synthetic UI data\n")

% Find entries that are of type char
isCharEntries = cellfun(@ischar,varargin);
indexOfCharEntries = find(isCharEntries);
charEntries = varargin(isCharEntries);

unknownEntries = ~ismember(charEntries,{'model','uLims','jLims','uErr','jErr','nSamp'});
if any(unknownEntries)
    unknownEntryCalls = string(charEntries(unknownEntries));
    error('createSyntheticUI:invalidData',"Invalid data! Unknown options: "...
        + join(unknownEntryCalls,", ")...
        + "\nOptions allowed: model, uLims, jLims, uErr, jErr, nSamp")
end

% Find locations of different parameter calls
modelCall = ismember(charEntries,{'model'});
jLimsCall  = ismember(charEntries,{'jLims'});
uLimsCall  = ismember(charEntries,{'uLims'});
jErrCall  = ismember(charEntries,{'jErr'});
uErrCall  = ismember(charEntries,{'uErr'});
nSampCall = ismember(charEntries,{'nSamp'});



%% Parse parameter calls
try
    if any(jLimsCall) % Current limits
        jLims = varargin{indexOfCharEntries(jLimsCall)+1};
    else
        jLims = [0.02 5];
    end
    
    if any(uLimsCall) % Voltage limits
        uLims = varargin{indexOfCharEntries(uLimsCall)+1};
    else
        uLims = [0 inf];
    end
    
    if any(jErrCall) % Standard deviation for current
        jErr = varargin{indexOfCharEntries(jErrCall)+1};
    else
        jErr = 0;
    end
    
    if any(uErrCall) % Standard deviation for voltage
        uErr = varargin{indexOfCharEntries(uErrCall)+1};
    else
        uErr = 0;
    end
    
    if any(nSampCall) % Standard deviation for voltage
        nSamp = varargin{indexOfCharEntries(nSampCall)+1};
    else
        nSamp = 200;
    end
    
catch ME
    if strcmp(ME.identifier,'MATLAB:badsubscript')
        error('Input must be given as name-value pairs')
    else
        rethrow(ME)
    end
end

%% Create electrolyzer model object

try
    if any(modelCall) 
        eModel = copy(varargin{indexOfCharEntries(modelCall)+1});
        
        if isempty(fieldnames(eModel.potentialFunc.Workspace.Variables)) % Are variables defined with the model?
            fprintf("\nNo variables provided with the model.\nDefault values used.\n")
            if strcmpi(eModel.type,'pem')
                Workspace = struct('Variables',struct('T',273.15+50,'pAn',2,'pCat',30));
            else % alkaline
                Workspace = struct('Variables',struct('T',273.15+50,'ps',2,'m',0.05));
            end
            eModel.setParams(Workspace);
        end
        
        if isempty(fieldnames(eModel.potentialFunc.Workspace.Coefficients)) % Are coefficients defined with the model?
            fprintf("\nNo fit coefficients provided with the model. \nDefault values used.\n")
            Workspace = struct('Coefficients',struct('alpha',0.5,'j0',1e-5,'r',0.1,'j_lim',1.5));
            eModel.setParams(Workspace);
        end
        
        if func.isEmpty(eModel.potentialFunc) % Are potential terms defined with the model?
            funcGiven = false;
            fprintf("\nModel provided but it did not have potential \nfunction defined. Default potentials used.\n")
            eModel.addPotentials('ocv','ohm','con');
        else
            funcGiven = true;
        end
        
    else
        fprintf("\nNo model provided. PEM electrolyzer and default values used for needed variables and coefficients.\n")
        funcGiven = false;
        eModel = electrolyzerModel();
        Variables = struct('T',273.15+50,'pAn',2,'pCat',30);
        Coefficients = struct('alpha',0.5,'j0',1e-5,'r',0.1,'j_lim',1.5);
        Workspace = struct('Variables',Variables,'Coefficients',Coefficients);
        eModel.setParams(Workspace);
        eModel.addPotentials('ocv','ohm','con');
    end
catch ME
    if strcmp(ME.identifier,'MATLAB:badsubscript')
        error('Input must be given as name-value pairs')
    else
        rethrow(ME)
    end
end

%% Create synthetic current and voltage vectors

% Get the wanted length for the vectors
if nargin >= 1 && isnumeric(varargin{1}) && isscalar(varargin{1})
    N = varargin{1};
else
    N = 20;
end

% If concentration overpotential is used, do not let current vector exceed
%   the limiting current value (j_lim)
if ismember('j_lim',fieldnames(eModel.potentialFunc.Workspace.Coefficients))
    jLims(2) = min(jLims(2),eModel.potentialFunc.Workspace.Coefficients.j_lim);
end

if funcGiven
    jmeas = linspace(jLims(1),jLims(2),10000);
    Umeas = eModel.calculate('current',jmeas);
else
    fprintf('\nActivation overpotential modelling properties:\n')
    fprintf('Activation voltage model: Buttler-Volmer equation\n')
    
    % Check that all the variables are scalars
    T = eModel.potentialFunc.Workspace.Variables.T;
    variableValues = struct2cell(eModel.potentialFunc.Workspace.Variables);
    for i = 1:length(variableValues)
        if numel(variableValues{i})>1
            error("Synthetic data with Buttler-Volmer equation can be created only with scalar variables.")
        end
    end
    
    % Calculate current and voltage vectors
    [F,R,n_e] = getConstants();
    f = R/(F*n_e);
    U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...50, Activation overpotential sweep limits
    
    Umeas = []; % "Measured" voltage
    jmeas = []; % "Measured" current based on Buttler-Volmer equation
    
    for ii = 1:length(U1)
        jmeastemp = eModel.potentialFunc.Workspace.Coefficients.j0*(exp(eModel.potentialFunc.Workspace.Coefficients.alpha/(f*T)*U1(ii))-exp((eModel.potentialFunc.Workspace.Coefficients.alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
        if jmeastemp >= jLims(2)
            break;
        elseif jmeastemp > jLims(1)
            U2 = eModel.calculate('current',jmeastemp); % Other potentials
            Umeastemp = U1(ii)+U2; % "Measured" voltage
            jmeas = [jmeas;jmeastemp];
            Umeas = [Umeas;Umeastemp];
        end
    end
end

% Take samples from the dense data vectors
Usamples = linspace(max(min(Umeas),uLims(1)),min(max(Umeas),uLims(2)),N)';
Udif = abs(Umeas-Usamples');
[~,ind] = min(Udif,[],1);
jmeassamp = jmeas(ind); % Final sampled current vector
Umeassamp = Umeas(ind); % Final sampled voltage vector



% Adding normal error with the given standard deviation to measurements
% The first value is kept errorless to avoid too radical change in current
% value that would hinder fitting performance.
jmeassamper = nan(length(jmeassamp),2);
Umeassamper = nan(length(jmeassamp),2);
for ii = 1:length(jmeassamp)
    jm = jmeassamp(ii).*(1 + randn(nSamp,1).*jErr);
    Um = Umeassamp(ii).*(1 + randn(nSamp,1).*uErr);
%     jm = jmeassamp(ii) + randn(nSamp,1).*jErr;
%     Um = Umeassamp(ii) + randn(nSamp,1).*uErr;
    jmeassamper(ii,:) = [mean(jm) std(jm)];
    Umeassamper(ii,:) = [mean(Um) std(Um)];
end

SynData = struct('voltage',Umeassamper,'current',jmeassamper);
FullData = struct('voltage',Umeas,'current',jmeas);

Workspace = eModel.potentialFunc.Workspace;

fprintf("\nSynthetic UI data creation finished.\n"...
    + "------------------------------------------------------------\n")
    

    

