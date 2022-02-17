classdef func < handle
    % FUNC Stores a structure-based function handle together with its workspace.
    %   
    %   The class FUNC enables simple combination of function handles
    %   maintaining their readability to the user. Objects of FUNC have a
    %   single function handle assigned to them even if some of the requred
    %   parameters contain data in vector form. FUNC manages this by having
    %   a separate Workspace structure that contains all the constants,
    %   variables and parameters in their respective structures. The
    %   function handle uses this Workspace as its only parameter and calls
    %   the required fields in the function handle. As there is only one
    %   input parameter, there is no requirement for the exact order of
    %   inputs for the function handle.
    %
    %   FUNC Properties:
    %       equation -- Function handle as an easily human-readable string.
    %       funcHandle -- The structure based function handle that uses
    %                       only Workspace as an input.
    %       Workspace -- The workspace structure containing substructures:
    %                       - Constants for system related constants, like
    %                           Faraday's constant or universal gas
    %                           constant.
    %                       - Variables for measured variable values, like
    %                           temperature or pressure that are needed for
    %                           the model.
    %                       - Parameters for system parameters that are
    %                           obtained from the parametrization, like
    %                           electric resistance or limiting current
    %                           density. Performing fit using
    %                           fitUI.m function automatically fills the
    %                           values for the parameters. Alternatively
    %                           they can be manually provided using
    %                           replaceInWorkspace method.
    %                       - Dependencies between any two variables. Every
    %                           time any Workspace variables are modified,
    %                           the dependencies are run to refresh the
    %                           dependent variables. Dependencies are
    %                           stored as strings of the commands and
    %                           function calls that are used to calculate
    %                           the dependent variables.
    %       Fitlims -- The structure for model parameter limits in
    %                       fitting. 1x3 cell array containing values
    %                       {low,start,high}
    %                       for each of the parameters. Any of the values
    %                       can be either numeric constants or MATLAB
    %                       equations in string form that output a numeric
    %                       constant using the independent variable of the
    %                       desired fit (marked with x).
    %
    %   FUNC Methods:
    %       calculate -- Calculates voltage from the UI curve based on
    %                       workspace and user-given values.
    %       copy -- Creates a copy of the object with a new handle.
    %       destructurize -- Destructurizes the structure-based function
    %                           handle to enable fitting with original
    %                           Matlab functions.
    %       removeFromWorkspace -- Remove values from the Workspace
    %                               structure. 
    %       replaceInWorkspace -- Replaces the values of existing
    %                               fields in Workspace structure.
    %       setFitlims -- Sets the values of Fitlims structure.
    %       setFuncHandle -- Sets the protected property of funcHandle
    %                           together with the property equation.
    %       setInWorkspace -- Sets new values to Workspace structure.
    %       viewWorkspace -- Outputs a human-readable report of the
    %                           contents of the Workspace structure.
    %       vectorify -- Destructurizes the structure-based function
    %                       handle to a form with vector inputs to enable
    %                       fitting with original Matlab functions.
    %   FUNC static methods:
    %       add -- Use addition to combine two func objects into a new one.
    %       createEmpty -- Creates an empty func object.
    %       isEmpty -- Determine if the func is empty.
    %       isWorkspace -- Determine if the given structure fulfills the
    %                       requirements for a Workspace.
    %
    %   See also FUNCTION_HANDLE, ELECTROLYZERMODEL
    
    properties (SetAccess = protected)
       equation; % Function handle in string form
       funcHandle; % Function handle of the function that uses only Workspace structure as an input
       Workspace; % Structure containing the Parameters, Variables, Constants and Dependencies between them:
%         Parameters; % Structure containing the parameters
%         Variables; % Structure containing the variables
%         Constants; % Structure containing the constants
%         Dependencies; % Structure containing the dependencies as strings
       Fitlims; % Limits for fitting the parameters in a structure containing fields [low start high]
    end
    
    %% Public methods
    methods
        function obj = func(funcHandle,Workspace,Fitlims)
            % FUNC Constructor method for the object
            %   Inputs:
            %       funcHandle -- The structure-based function handle.
            %       Workspace -- The workspace structure.
            %       Fitlims (optional) -- The fit limit structure
            obj.setFuncHandle(funcHandle);
            if func.isWorkspace(Workspace)
                obj.Workspace = Workspace;
            else
                error("Workspace structure of object func should contain exclusively fields 'Parameters', 'Variables', 'Constants' or 'Dependencies'.")
            end
            
            obj.Fitlims = struct();
            if exist("Fitlims","var")==1&&~isempty(Fitlims)&&isstruct(Fitlims) % Fit limits provided
                obj.setFitlims(Fitlims);
            end    
            
            obj.refreshWorkspace;
        end
        

        function setFuncHandle(obj,newFuncHandle)
            % SETFUNCHANDLE Sets the funcHandle and equation properties.
            obj.funcHandle = newFuncHandle;
            obj.equation = obj.getEquation;
        end
        
        
        function setInWorkspace(obj,SetStruct)
            % SETINWORKSPACE  Sets new values in the Workspace structure.
            %  Parameters should be provided as a Workspace-compatible 
            %  structure. The method prioritizes given values over the
            %  existing ones in case of conflict.
            %
            %  See also FUNC.ISWORKSPACE, MERGESTRUCTS
            if func.isWorkspace(SetStruct)
                OldWorkspace = obj.Workspace;
                obj.Workspace = mergeStructs(obj.Workspace,SetStruct);
                
                obj.refreshWorkspace(OldWorkspace);
            else
                error('Values to be set have to be given as a single Workspace structure')
            end
            
        end
        
        
        function replaceInWorkspace(obj,varargin)
            % REPLACEINWORKSPACE  Replaces values in the Workspace structure.
            %  Values should be provided as name-value pairs or as a
            %  structure (either Workspace compatible or not). Replaces
            %  values only for existing fields with the same name as given.
            %  The method acts recursively inside Workspace.
            %  
            %  Note: If a dependency is defined for a parameter to be
            %  replaced, the dependency is erased, unless option 'rebuild'
            %  is used. Rebuild option is used when REPLACEINWORKSPACE
            %  is used to input fit results when introducing values to a
            %  copied func object. The benefits of the option are in
            %  situations when there are dependencies defined that do not
            %  directly alter the parameter they are related to, but
            %  execute some other code, like a warning. If 'rebuild' is
            %  used when there are dependencies affecting a parameter that
            %  is supposed to be changed manually, the dependency
            %  overwrites any changes made.
            %
            % See also ADDVALUESTOSTRUCT
            
            rebuildCall = strcmpi(varargin,'rebuild');
            varargin(rebuildCall) = [];
            
            if isempty(varargin{1})
                return;
            elseif length(varargin) == 1 && isstruct(varargin{1}) % Single structure input
                Struct = varargin{1};
                if func.isWorkspace(Struct) % Workspace structure input
                    fields = fieldnames(Struct);
                    for i = 1:length(fields)
                        if any(rebuildCall)
                            replaceInWorkspace(obj,Struct.(fields{i}),'rebuild')
                        else
                            replaceInWorkspace(obj,Struct.(fields{i}))
                        end
                    end
                    return;
                else % Non-workspace structure input
                    valuesToReplace = Struct;
                end
            elseif ~mod(length(varargin),2) % Input as name-value pairs
                for i = 1:2:length(varargin)
                    valuesToReplace.(varargin{i}) = varargin{i+1};
                end
            else
                error('Parameters to be replaced have to be either set as a single structure or as name-value pairs')
            end
            
            % Means of avoiding overwriting existing dependencies
            if ismember('Dependencies',fieldnames(obj.Workspace))
                paramNames = fieldnames(valuesToReplace);
                dependentParamNames = paramNames(isfield(obj.Workspace.Dependencies,paramNames));
                tempDependencies = obj.Workspace.Dependencies; % Store the 
                    % dependencies from the old Workspace to be manually
                    % added to the new Workspace after replacing the
                    % parameters.
                if ~any(rebuildCall) % Erase dependencies for the replaced parameters.
                    % Enables manual input of dependent variables so that
                    % they are no longer dependent.
                    tempDependencies = rmfield(tempDependencies,dependentParamNames);
                end
            else
                tempDependencies = [];
            end
            
            OldWorkspace = obj.Workspace;
            obj.Workspace = addValuesToStruct(obj.Workspace,valuesToReplace);
            if ~isempty(tempDependencies)
                obj.Workspace.Dependencies = tempDependencies;
            end
            
            obj.refreshWorkspace(OldWorkspace);
        end
        
        function removeFromWorkspace(obj,varargin)
            % REMOVEFROMWORKSPACE  Removes parameters from the Workspace structure.
            %  Input has to be provided either as a list of
            %  strings/character vectors or as a cell array of
            %  strings/character vectors containing the names of the
            %  parameters to be removed. The method acts recursively inside
            %  the Workspace.
            
            if length(varargin) == 1 && iscell(varargin{1})
                valuesToRemove = varargin{1};
            else
                valuesToRemove = varargin;
            end
            
            OldWorkspace = obj.Workspace;
            
            for i = 1:length(valuesToRemove)
                if ismember(valuesToRemove{i},fieldnames(obj.Workspace.Constants))
                    obj.Workspace.Constants = rmfield(obj.Workspace.Constants,valuesToRemove{i});
                end
                if ismember(valuesToRemove{i},fieldnames(obj.Workspace.Parameters))
                    obj.Workspace.Parameters = rmfield(obj.Workspace.Parameters,valuesToRemove{i});
                end
                if ismember(valuesToRemove{i},fieldnames(obj.Workspace.Variables))
                    obj.Workspace.Variables = rmfield(obj.Workspace.Variables,valuesToRemove{i});
                end
                if ismember(valuesToRemove{i},fieldnames(obj.Workspace.Dependencies))
                    obj.Workspace.Dependencies = rmfield(obj.Workspace.Dependencies,valuesToRemove{i});
                end
            end
            
            obj.refreshWorkspace(OldWorkspace);
        end
        
        function setFitlims(obj,varargin)
            % SETFITLIMS Sets the property fitLims.
            %   Set fitting limit values for the model parameters.
            %   Fitlims structure contains one field for each parameter,
            %   named identical to the parameter. Each field contains a
            %   1x3 cell array with values {low,start,high} for lower
            %   limit, starting point and higher limit, respectively.
            %
            %   Recognized input is a Structure with the right format or
            %   name value pairs with the parameter name followed by the
            %   limit matrix.
            %
            %   Limits can be either as:
            %       - numeric scalars or
            %       - MATLAB equations in string form each outputing a
            %       numeric scalar based on the dependent variable values.
            %       Dependent variable has to be signed with 'x' in the
            %       equations.
            %
            %   Example of the Fitlim structure:
            %       struct('a',{{1,2,3}},... % Numeric scalars
            %              'b',{{'min(x)','avg(x)','max(x)'}}) % Equations
            
            if isempty(varargin{1})
                return;
            elseif length(varargin) == 1 && isstruct(varargin{1})
                Struct = varargin{1};
            elseif mod(nargin,2)
                for i = 1:2:length(varargin)
                    Struct.(varargin{i}) = varargin{i+1};
                end
            else
                error('Fit limits have to be set either as a single structure or as name-value pairs')
            end
            
            fn = fieldnames(Struct);
            for i = 1:numel(fn)
                if ~isequal(size(Struct.(fn{i})),[1 3])||~iscell(Struct.(fn{i}))
                    error("Given Fitlimit structure has invalid values. Provide the limits as a 1x3 cell array {low,start,high}")
                end
            end
            obj.Fitlims = mergeStructs(obj.Fitlims,Struct,'warn_duplicates',false);
        end

        
        
        function result = calculate(obj,varargin)
            % CALCULATE Calculates the result from the function handle.
            %   Variables input as name-value pairs are used for the
            %   calculation and any variable required by the function that is
            %   not input to CALCULATE is looked for from the Workspace of the
            %   object.
            %
            %   Examples of usage:
            %
            %   results = obj.CALCULATE calculates voltage values based on
            %               the Workspace structure. Also current has to be
            %               included in the workspace
            %
            %   [] = obj.CALCULATE('current',j) calculates voltage values
            %               based on the current values provided in vector
            %               j. Other needed variables are taken from the
            %               Workspace.
            %
            %   [] = obj.CALCULATE('var1',var1,'var2',var2) calculates
            %               voltage values based on the provided values and
            %               obtains the rest of the needed variables from
            %               the Workspace.
            %
            %   [] = obj.CALCULATE('Workspace',Workspace) calculates
            %               voltage values based on the provided
            %               Workspace-structure only.
            %
            %   [] = obj.CALCULATE(__,'detectChanges',b) enables user to
            %               decide whether dependencies that are related
            %               to the change in some Workspace parameter are
            %               evaluated or not. This kind of dependencies are
            %               mostly warnings about unaccounted temperature
            %               dependency in some models. Options for b:
            %                   true -- Detect changes and evaluate change
            %                       dependent dependencies. This may slow
            %                       down calculation in the case of looped
            %                       calculation
            %                   false -- Do not evaluate change dependent
            %                       dependencies. 
            %
            %   Note: Providing full Workspace or individual parameters to
            %   CALCULATE method does not modify the Workspace stored in
            %   the func object.
            
            
            if ~mod(nargin,2)
                error("Variables have to be given as name-value pairs")
            end
            
            % Check if change detection was specified by the user
            detectChangesCall = strcmpi(varargin,'detectChanges');
            if any(detectChangesCall)
                index = [false detectChangesCall(1:end-1)]; % index after the call
                detectChanges = varargin{index};
                varargin = varargin(~(detectChangesCall & index));
            else % use default
                detectChanges = true;
            end
            
            
            if isequal(obj.funcHandle,@pass)
                error("No function defined. Results cannot be calculated.")
            end
            
            % Check if Workspace was given as an input
            workspaceCall = strcmpi(varargin,'Workspace');
            if any(workspaceCall)
                % Create temporary Workspace based only on the given
                % Workspace
                index = [false workspaceCall(1:end-1)]; % index after the call
                TempWorkspace = func.createTempWorkspace(varargin{index});
            else
                % Create a temporary workspace that includes the user input
                % variables in addition to the ones found from the func 
                % Workspace.                
                TempWorkspace = func.createTempWorkspace(obj.Workspace,varargin);
            end
            
            % Refresh dependencies
            if detectChanges
                % Refresh dependencies that relate to change state of some
                % Workspace parameters (mostly warnings)
                TempWorkspaceRefreshed = func.refresh(TempWorkspace,obj.Workspace);
            else
                % Refresh only dependencies that are directly related to
                % the values of Workspace parameters
                TempWorkspaceRefreshed = func.refresh(TempWorkspace);
            end
            
            % The result is calculated if the TempWorkspace contains all
            % necessary values for the calculation.
            try
                result = obj.funcHandle(TempWorkspaceRefreshed);
            catch ME
                if strcmp(ME.identifier,'MATLAB:nonExistentField')
                    msg = "Result cannot be calculated because one or " +...
                            "more necessary workspace entries are missing. " +...
                            "\n\nError message: "+...
                            string(ME.message) +...
                            "\n\nEither add the missing entries to the " +...
                            "workspace of the func object or input them " +...
                            "to the calculate method as name-value pairs.";
                    newME = MException('MATLAB:func:incompleteWorkspace',compose(msg));
                    throw(newME)
                else
                    rethrow(ME)
                end
            end
        end
        
        
        function [destructurizedFuncHandle,...
                parametersNamesForFit,...
                problemVariableNames,...
                problemVariables] = destructurize(obj,...
                                                  independentVariableName)
            % DESTRUCTURIZE Destructurizes the structure-based function handle.
            %   Modifies the structure-based function handle to a function
            %   handle with all the constants and scalar variables
            %   calculated in, and the vector variables (problem variables)
            %   and all the fit parameters used as separate input
            %   parameters.
            %   Destructurizing enables usage of original Matlab functions
            %   like fit and particleswarm.
            %
            %   The name of the independent variable has to be provided for
            %   the correct ordering of the variables. 
            %
            %   Inputs:
            %       independentVariableName -- Name of the independent
            %                                   variable that is kept as
            %                                   the last parameter to the
            %                                   destructurized funtion
            %                                   handle.
            %   Outputs:
            %       destructurizedFuncHandle -- The destructurized function
            %                                   handle.
            %       parameterNamesForFit -- Lists of the undefined 
            %                                   parameter names as a cell 
            %                                   array.
            %       problemVariableNames -- List of the names of variables
            %                               with vector data as a cell 
            %                               array.
            %       problemVariables -- Cell array of the data of the
            %                           variables listed in method
            %                           problemVariableNames.
            
            % Constants
            if ~isCompleteStruct(obj.Workspace.Constants)
                error('Missing constant values. Not able to destructurize function handle for fitting.')
            end
            %
            % Parameters
            parametersNamesForFit = {};
            constantParameterNames = {};
            constantParameters = [];
            if ~isempty(obj.Workspace.Parameters)
                allParameterNames = fieldnames(obj.Workspace.Parameters);
                for i = 1:length(allParameterNames)
                    parameter = obj.Workspace.Parameters.(allParameterNames{i});
                    if contains(obj.getEquationBody,['Workspace.Parameters.' allParameterNames{i}])
                        % If the parameter is used in the equation
                        parametersNamesForFit = [parametersNamesForFit;...
                                                   allParameterNames{i}];
                    else % If the parameter is not used in the equation
                        if ~isempty(parameter)
                            constantParameterNames = [constantParameterNames;...
                                                        allParameterNames{i}];
                            constantParameters = [constantParameters;...
                                                    parameter];
                        else
                            warning("There is a parameter whose name is not found from the equation but no value has been assigned to it.")
                        end
                    end
                end
            else
                warning("No parameters defined for the function.")
            end

            % Variables
            allVariableNames = fieldnames(obj.Workspace.Variables);
            if ismember(independentVariableName,allVariableNames)
                independentVariable = obj.Workspace.Variables.(independentVariableName);
            else
                error('Inependent variable is not included in the workspace. Not able to destructurize function handle for fitting.')
            end
            problemVariableNames = {};
            problemVariables = {};
            nonProblemVariableNames = {};
            nonProblemVariables = [];
            for i = 1:length(allVariableNames)
                variableName = allVariableNames{i};
                variable = obj.Workspace.Variables.(variableName);
                if isstring(variable)||ischar(variable)
                    continue
                elseif strcmp(allVariableNames{i},independentVariableName)
                    continue
                elseif isempty(variable)
                    error('One or more variables missing. Not able to destructurize function handle for fitting.')
                elseif length(variable) == 1
                    nonProblemVariableNames = [nonProblemVariableNames;...
                                               variableName];
                    nonProblemVariables = [nonProblemVariables;variable];
                else
                    problemVariableNames = [problemVariableNames;...
                                            variableName];
                    problemVariables = [problemVariables;variable];
                end                    
            end

            % Create destructurized function handle

            % Constants
            Constants = obj.Workspace.Constants;

            % Parameters
            Parameters = struct();
            for i = 1:length(parametersNamesForFit)
                Parameters.(parametersNamesForFit{i}) = sym(parametersNamesForFit{i});
            end
            %
            for i = 1:length(constantParameterNames)
                Parameters.(constantParameterNames{i}) = constantParameters(i);
            end
           
            % Variables
            for i = 1:length(nonProblemVariableNames)
                Variables.(nonProblemVariableNames{i}) = nonProblemVariables(i);
            end
            %
            for i = 1:length(problemVariableNames)
                Variables.(problemVariableNames{i}) = sym(problemVariableNames{i});
            end
            
            Variables.(independentVariableName) = sym(independentVariableName);
            
            % Build new function handle in symbolic form
            symbolicFunction = obj.funcHandle(struct('Constants',Constants,...
                                                     'Variables',Variables,...
                                                     'Parameters',Parameters));
            
            % Modify back to function handle
            varsList = [parametersNamesForFit;...
                        problemVariableNames;...
                        independentVariableName];
            destructurizedFuncHandle = matlabFunction(symbolicFunction,...
                'Vars',varsList);
        end
        
        
        
        
        function [vectorifiedFuncHandle,...
                parametersNamesForFit,...
                problemVariableNames,...
                problemVariables] = vectorify(obj,...
                                              independentVariableName)
            % VECTORIFY returns a function handle with one vector for the
            %   parameters transformed from input function handle with separate
            %   parameters that are organized in order
            %   (parameters, problem variables, independent variable).
            %   Inputs:
            %       independentVariableName -- Name of the independent
            %                                   variable that is kept as
            %                                   the last parameter to the
            %                                   destructurized funtion
            %                                   handle.
            %   Outputs:
            %       vectorifiedFuncHandle -- Function handle modified to
            %                                   take cell array inputs for
            %                                   parameters instead of
            %                                   listing them separately.
            %       parameterNamesForFit -- Lists of the undefined 
            %                                   parameter names as a cell 
            %                                   array.
            %       problemVariableNames -- List of the names of variables
            %                               with vector data as a cell 
            %                               array.
            %       problemVariables -- Cell array of the data of the
            %                           variables listed in method
            %                           problemVariableNames.
            %
            %   See also FUNC.DESTRUCTURIZE
            
            [destructurizedFuncHandle,...
                parametersNamesForFit,...
                problemVariableNames,...
                problemVariables] = destructurize(obj,...
                                                  independentVariableName);
            
            nParams = numel(parametersNamesForFit);
            nProbVars = numel(problemVariableNames);
                                              
            % Creating symbolic variables for the function handle
            x = num2cell(sym('x', [1 nParams])); % fit parameters
            y = num2cell(sym('y', [1 nProbVars])); % problem variables
            syms independentVariable; % Independent variable input
            
            % Calculating the symbolic result
            z = destructurizedFuncHandle(x{:},y{:},independentVariable);
            
            % Creating a modified function handle from the symbolic result
            vectorifiedFuncHandle = matlabFunction(z,'Vars',[{cell2sym(x)},y(:)',{independentVariable}]);
            
        end
        
        
        
        
        function childFunc = copy(obj,varargin)
            % COPY Creates an independent copy of the object.
            %   Changing properties of the copy doesn't affect the
            %   original.
            %
            % With input parameter 'empty' creates a copy that has no
            %  funcHandle but keeps the Workspace content.
            
            if nargin == 1
                for i = 1:height(obj)
                    for j = 1:width(obj)
                        childFunc(i,j) = func(obj(i,j).funcHandle,obj(i,j).Workspace,obj(i,j).Fitlims);
                    end
                end
            elseif strcmp(varargin{1},'empty')
                for i = 1:numel(obj)
                    for j = 1:width(obj)
                        childFunc(i,j) = func(@pass,obj(i,j).Workspace,obj(i,j).Fitlims);
                    end
                end
            else
                error('Unknown input')
            end
        end
        
        
        function varargout = viewWorkspace(obj)
            % VIEWWORKSPACE Outputs a tabular report of the Workspace.
            
            fieldNameConst = fieldnames(obj.Workspace.Constants);
            fieldNameVars = fieldnames(obj.Workspace.Variables);
            fieldNameParam = fieldnames(obj.Workspace.Parameters);
            fields = [strcat(cell(size(fieldNameConst)),'Constants'),fieldNameConst;...
                         strcat(cell(size(fieldNameVars)),'Variables'),fieldNameVars;...
                         strcat(cell(size(fieldNameParam)),'Parameters'),fieldNameParam];
            fieldName = fields(:,2);
            description = cell(size(fieldName));
            valueMin = nan(size(fieldName));
            valueMean = nan(size(fieldName));
            valueMax = nan(size(fieldName));
            standardDeviation = nan(size(fieldName));
            for i = 1:length(fieldName)
                value = obj.Workspace.(fields{i,1}).(fields{i,2});
                if strcmpi(fields(i,2),'electrolyte')
                    descriptionText = ": numeric helper variable for alkaline";
                    valueMean(i) = value(:,1);
                elseif ismember('Dependencies',fieldnames(obj.Workspace)) && any(ismember(fieldName(i),fieldnames(obj.Workspace.Dependencies)))
                    descriptionText = ": dependent";
                    if isempty(value)
                        descriptionText = strcat(descriptionText,", no value assigned");
                    elseif length(value(:,1)) == 1
                        valueMean(i) = value(:,1);
                        descriptionText = strcat(descriptionText," scalar");
                    else
                        valueMin(i) = min(value(:,1));
                        valueMax(i) = max(value(:,1));
                        valueMean(i) = mean(value(:,1));
                        descriptionText = strcat(descriptionText,strcat(" vector of length ",num2str(length(value(:,1)))));
                    end
                elseif ~isempty(value)
                    if length(value(:,1)) == 1
                        valueMean(i) = value(:,1);
                        descriptionText = ": scalar";
                        if length(value(1,:)) == 2
                            standardDeviation(i) = value(1,2);
                            descriptionText = strcat(descriptionText," with confidence bounds");
                        else
                            descriptionText = strcat(descriptionText," without confidence bounds");
                        end
                    else
                        valueMin(i) = min(value(:,1));
                        valueMax(i) = max(value(:,1));
                        valueMean(i) = mean(value(:,1));
                        descriptionText = strcat(": vector of length ",num2str(length(value(:,1))));
                        if ~kstest(value(:,1)) % If data is normally distributed
                            standardDeviation(i) = std(value(:,1));
                            descriptionText = strcat(descriptionText,", normally distributed");
                        elseif length(value(1,:)) == 2
                            descriptionText = strcat(descriptionText," with individual values of standard deviation");
                        end
                    end
                else % No value assigned
                    descriptionText = ": no values assigned";    
                end
                description{i} = strcat(fields{i,1}(1:end-1),descriptionText);
            end
            report = table(description,valueMean,valueMin,valueMax,standardDeviation);
            report.Properties.RowNames = fieldName;
            
            if nargout >= 1
                varargout{1} = report;
            else
                disp(report)
            end
        end
        
    end
    
    
    %% Public static methods
    
    methods(Static, Access = public)
        function newFunc = add(func1,func2)
            % ADD Add two func objects together.
            %   Combining the workspaces and uses addition for combining
            %   their function handles.
            %
            %   See also FUNC, MERGESTRUCTS
            
            if ~isa(func1,'func')||~isa(func2,'func')
                error("Functions to be added should both be of type 'func'")
            end
            funcs = {func1,func2};
            emptyFuncs = func.isEmpty(funcs);
            if any(emptyFuncs)
                funcToAdd = funcs(~emptyFuncs);
                newFuncEquation = ['@(Workspace) ' funcToAdd{:}.getEquationBody];
            else
                newFuncEquation = ['@(Workspace) ' func1.getEquationBody ' + ' func2.getEquationBody];
            end
            newFuncHandle = str2func(newFuncEquation);
            
            NewWorkspace = mergeStructs(func1.Workspace,func2.Workspace);
            
            NewFitlims = mergeStructs(func1.Fitlims,func2.Fitlims);
            
            newFunc = func(newFuncHandle,NewWorkspace,NewFitlims);
            
%             newFunc.refreshWorkspace;
            
        end
        
        function emptyFunc = createEmpty
            % CREATEEMPTY Create an empty func object.
            %  Creates a func object with a function handle that does
            %  nothing, and an empty, but workspace-compatible Workspace
            %  structure.
            emptyFunc = func(@pass,...
                struct('Constants',struct(),...
                'Variables',struct(),...
                'Parameters',struct(),...
                'Dependencies',struct()));
        end
        
        function b = isEmpty(varargin)
            % ISEMPTY Evaluates if the given func objects are empty
            if length(varargin) == 1 && iscell(varargin{1})
                funcs = varargin{1};
            else
                funcs = varargin;
            end
            
            b = false(size(funcs));
            
            for i = 1:length(funcs)
                funcToEval = funcs{i};
                if strcmp(funcToEval.equation,'pass')
                    b(i) = true;
                end
            end
        end
        
        function b = isWorkspace(Struct)
            % ISWORKSPACE Evaluates Workspace-compatibility of a structure.            
            b = all(ismember(fieldnames(Struct),...
                             {'Parameters';...
                             'Variables';...
                             'Constants';...
                             'Dependencies'}));
            for fn = transpose(fieldnames(Struct))
                b = b && isstruct(Struct.(fn{:}));
            end
        end
    end
    
    %% Private dynamic methods
    methods (Access = private)
        
        function equationStr = getEquation(obj)
            % GETEQUATION Creates human-readable equation string.
            %   Erases all occurences of Workspace from the equation body
            %   leaving a string with human-readable equation structure.
            equationStr = erase(obj.getEquationBody,...
                {'Workspace.Parameters.',...
                'Workspace.Variables.',...
                'Workspace.Constants.'});
        end
        
        function equationBody = getEquationBody(obj)
            % GETEQUATIONBODY Outputs only the equation body.
            %  Removes the @(Workspace) from the function handle.
            equationBody = erase(func2str(obj.funcHandle),'@(Workspace)');
        end
        
        function refreshWorkspace(obj,varargin)
            % REFRESHWORKSPACE Recalculates dependent Workspace parameters.
            if ~isempty(varargin)
                obj.Workspace = func.refresh(obj.Workspace,varargin{1});
            else
                obj.Workspace = func.refresh(obj.Workspace);
            end
            
            
            % Determine if there are undefined parameters without fit limits
            % First find all the parameters that don't have an value
            fnParams = fieldnames(obj.Workspace.Parameters);
            emptyFields = false(size(fnParams));
            for fni = 1:numel(fnParams)
                emptyFields(fni) = isempty(obj.Workspace.Parameters.(fnParams{fni}));
            end
            fnEmpty = fnParams(emptyFields); % Field names of empty parameters
            
            if any(~ismember(fnEmpty,fieldnames(obj.Fitlims)))
                warningmsg = "Workspace contains parameters with no"...
                    + " set fit limits. Consider setting the limits with"...
                    + " the setFitlims method before using this func"...
                    + " object for fitting.";
                warning(warningmsg)
            end
        end
            
        
    end
    

    %% Private static methods
    
    methods (Access = private, Static)
        function TempWorkspace = createTempWorkspace(Struct,replacingArguments)
            % CREATETEMPWORKSPACE Creates a copy of the Workspace structure
            %   but without values for standard deviation. Replaces the
            %   values for given fields (replacingArguments, cell array of
            %   name-value pairs)
            
            fn = fieldnames(Struct);
            
            if isempty(replacingArguments)
                replArgNames = {};
            else
                replArgNames = replacingArguments(1:2:end);
                replArgs = replacingArguments(2:2:end);
            end
            
            TempWorkspace = Struct;
            
            for i = 1:length(fn)
                if strcmp(fn{i},'Dependencies')
                    TempWorkspace.(fn{i}) = Struct.(fn{i});
                elseif isstruct(Struct.(fn{i}))
                    TempWorkspace.(fn{i}) = func.createTempWorkspace(Struct.(fn{i}),replacingArguments);
                elseif isempty(Struct.(fn{i}))
                    replIndex = ismember(replArgNames,fn{i});
                    if any(replIndex)
                        TempWorkspace.(fn{i}) = replArgs{replIndex};
                    else
                        continue;
                    end
                elseif length(Struct.(fn{i})(1,:)) == 1
                    TempWorkspace.(fn{i}) = Struct.(fn{i});
                else
                    fieldSize = size(Struct.(fn{i}));
                    if fieldSize(2)>2
                        warningMsg = strcat("It seems that the Workspace variable "...
                            ,string(fn{i})," has more than two columns (",...
                            string(fieldSize(1)),"x",string(fieldSize(2)),...
                            "). Consider providing the data as the first column of the variable matrix. ",...
                            "The second column should be preserved for standard deviation, if applicable. ",...
                            "The method func.calculation takes into acount only the first column of the vector and regards that as the given dataset.");
                        warning(warningMsg)
                    end
                    TempWorkspace.(fn{i}) = Struct.(fn{i})(:,1);
                end
            end
        end
        
        function refreshedWorkspace = refresh(Workspace,varargin)
            % REFRESH Recalculates dependent Workspace parameters.
            if nargin >= 2
                OldWorkspace = varargin{1};
                changed = func.detectChanges(Workspace,OldWorkspace);
            else
                changed = func.detectChanges(Workspace,Workspace);
            end
            
            if ismember('Dependencies',fieldnames(Workspace))
                fn = fieldnames(Workspace.Dependencies);
                for i = 1:numel(fn)
                    try
                        warning off backtrace
                        eval(Workspace.Dependencies.(fn{i}))
                        warning on backtrace
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:InputParser:ArgumentFailedValidation')||strcmp(ME.identifier,'MATLAB:undefinedVarOrClass')
                            % Some variable needed in a function call in
                            % dependencies is not defined. This causes the
                            % function call to throw an error, which is
                            % unnecessary for this object
                            continue;
                        elseif strcmp(ME.identifier,'MATLAB:nonExistentField')
                            warningmsg = "Some dependencies could not be " + ...
                                "evaluated because of missing variables. " + ...
                                "Ensure that all the needed variables are " + ...
                                "defined in the Workspace before using " + ...
                                "this func object further.\n\n" + ...
                                ME.message;
                            warning('refreshWorkspace:nonExistentField',warningmsg)
                        else
                            rethrow(ME)
                        end
                    end
                end
            end
            refreshedWorkspace = Workspace;
        end
        
        
        function changed = detectChanges(NewWorkspace,OldWorkspace)
            % DETECTCHANGES Detects differences between two Workspaces
            if func.isWorkspace(NewWorkspace) && func.isWorkspace(OldWorkspace) 
                [common, dNew, ~] = comp_struct(NewWorkspace,OldWorkspace);
                compResults = {common,dNew};
                changed = struct();
                for i = 1:2
                    S = compResults{i};
                    if isstruct(S)
                        for fn1 = fieldnames(S)'
                            for fn2 = fieldnames(S.(fn1{:}))'
                                changed.(fn2{:}) = (i~=1);
                            end
                        end
                    end
                end      
            else
                error("One or both of the given structures are not Workspace-compatible")
            end
        end
        
    end
    
    
end