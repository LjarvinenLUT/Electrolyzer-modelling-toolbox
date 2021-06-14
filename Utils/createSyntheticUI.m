function [SynData,FullData,Workspace] = createSyntheticUI(varargin)
% CREATESYNTHETICUI Creates synthetic Voltage-Current pairs from an
%   electrolyzer UI curve.
%
%   [synData,Workspace] = CREATESYNTHETICUI() creates synthetic UI curve
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
%           - 'func', enables user-defined func object to be used for
%               creation of the synthetic UI curve.
%           - 'workspace', enables user to define all the variable and
%               coefficient values required for the creation of the
%               synthetic UI curve.
%           - 'jLims', enables defining of current limits for the synthetic
%               UI data. Input should be in a form of 1x2 vector 
%               [iLow iHigh].
%           - 'jSigma', enables defining of standard deviation for current
%               measurements.
%           - 'uLims', enables defining of voltage limits for the synthetic
%               UI data. Input should be in a form of 1x2 vector 
%               [uLow uHigh].
%           - 'uSigma', enables defining of standard deviation for voltage
%               measurements.
%           - 'type', electrolyzer type. Possible options: "PEM" and
%               "alkaline". Default: "PEM"
%           - 'electrolyte', electrolyte type for alkaline electrolyzers.
%               Possible options: "KOH" and "NaOH". Default: "KOH"
%
%   Output:
%       - SynData -- Structure containing sampled data in fields 'current' 
%           and 'voltage', both of them being Nx2 matrices:
%           (first column = data, second column = standard deviation).
%       - FullData -- Structure containing full data in fields 'current' 
%           and 'voltage', both of them being Nx2 matrices:
%           (first column = data, second column = standard deviation).
%       - Workspace -- Structure including the variables, constants and
%           coefficients used for creating the synthetic data.
%
%   See also FITUI, FUNC, ELECTROLYZERMODEL

fprintf("\nCreating synthetic UI data\n")

% Find entries that are of type char
isCharEntries = cellfun(@ischar,varargin);
indexOfCharEntries = find(isCharEntries);
charEntries = varargin(isCharEntries);

% Find locations of different parameter calls
funcCall = ismember(charEntries,{'func'});
workspaceCall = ismember(charEntries,{'workspace'});
jLimsCall  = ismember(charEntries,{'jLims'});
uLimsCall  = ismember(charEntries,{'uLims'});
jSigmaCall  = ismember(charEntries,{'jSigma'});
uSigmaCall  = ismember(charEntries,{'uSigma'});
typeCall = ismember(charEntries,{'type'});
electrolyteCall = ismember(charEntries,{'electrolyte'});

%% Parse parameter calls
try
    if any(typeCall) % Electrolyzer type
        type = varargin{indexOfCharEntries(typeCall)+1};
    else
        type = "pem";
        disp("No electrolyzer type provided; Synthetic UI created for PEM.")
    end
    
    if any(electrolyteCall) % Electrolyte for alkaline electrolysis
        electrolyte = varargin{indexOfCharEntries(electrolyteCall)+1};
    else
        electrolyte = "KOH";
        if strcmp(type,'alkaline')
            disp("Electrolyte not specified for alkaline electrolyzer; KOH used as a default.")
        end
    end
    
    if any(jLimsCall) % Current limits
        jLims = varargin{indexOfCharEntries(jLimsCall)+1};
    else
        jLims = [0 inf];
    end
    
    if any(uLimsCall) % Voltage limits
        uLims = varargin{indexOfCharEntries(uLimsCall)+1};
    else
        uLims = [0 inf];
    end
    
    if any(jSigmaCall) % Standard deviation for current
        jSigma = varargin{indexOfCharEntries(jSigmaCall)+1};
    else
        jSigma = 0;
    end
    
    if any(uSigmaCall) % Standard deviation for voltage
        uSigma = varargin{indexOfCharEntries(uSigmaCall)+1};
    else
        uSigma = 0;
    end
    
catch ME
    if strcmp(ME.identifier,'MATLAB:badsubscript')
        error('Input must be given as name-value pairs')
    else
        rethrow(ME)
    end
end

%% Create electrolyzer model object
eModel = electrolyzerModel('type',type,'electrolyte',electrolyte);

try
    if any(funcCall) 
        funcGiven = true;
        eModel.addPotentials(varargin{indexOfCharEntries(funcCall)+1});
        if any(workspaceCall)
            eModel.setParams(varargin{indexOfCharEntries(workspaceCall)+1});
        end
    elseif any(workspaceCall)
        funcGiven = false;
        Workspace = varargin{indexOfCharEntries(workspaceCall)+1};
        eModel.setParams(Workspace);
        coeffNames = fieldnames(Workspace.Coefficients);
        potentialNames = {'ocv','ohm','con'};
        includedPotentials = [true ismember({'r','j_lim'},coeffNames)];
        eModel.addPotentials(potentialNames{includedPotentials});
    else
        funcGiven = false;
        if strcmp(type,"pem") % PEM
            Variables = struct('T',273.15+50,'pAn',2,'pCat',30);
        else % Alkaline
            Variables = struct('T',273.15+50,'p',2,'m',7.64);
        end
        Workspace = struct('Variables',Variables,'Coefficients',struct('alpha',0.5,'j0',1e-5,'r',0.1,'j_lim',1.5));
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
    fprintf('\nActivation overpotential calculation properties:\n')
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
iii = 1;
for ii = 1:N
    Udif = abs(Umeas-Usamples(ii));
    [~,ind] = min(Udif);
    if iii == 1 || (iii > 1 && abs(jmeas(ind)-jmeassamp(iii-1))>0.01)
        jmeassamp(iii,1) = jmeas(ind); % Final sampled current vector
        Umeassamp(iii,1) = Umeas(ind); % Final sampled voltage vector
        iii = iii+1;
    end
end

% Adding normal error with the given standard deviation to measurements
% The first value is kept errorless to avoid too radical change in current
% value that would hinder fitting performance.
jmeassamper = [jmeassamp(1) 0;nan(length(jmeassamp)-1,2)];
Umeassamper = [Umeassamp(1) 0;nan(length(jmeassamp)-1,2)];
for ii = 2:length(jmeassamp)
%     jm = jmeassamp(ii).*(1 + randn(200,1).*jSigma);
%     Um = Umeassamp(ii).*(1 + randn(200,1).*uSigma);
    jm = jmeassamp(ii) + randn(200,1).*jSigma;
    Um = Umeassamp(ii) + randn(200,1).*uSigma;
    jmeassamper(ii,:) = [mean(jm) std(jm)];
    Umeassamper(ii,:) = [mean(Um) std(Um)];
end

SynData = struct('voltage',Umeassamper,'current',jmeassamper);
FullData = struct('voltage',Umeas,'current',jmeas);
    

    

