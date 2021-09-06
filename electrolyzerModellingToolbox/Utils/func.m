classdef func < handle
    % FUNC Stores a structure-based function handle together with its workspace.
    %   
    %   The class FUNC enables simple combination of function handles
    %   maintaining their readability to the user. Objects of FUNC have a
    %   single function handle assigned to them even if some of the requred
    %   parameters contain data in vector form. FUNC manages this by having
    %   a separate Workspace structure that contains all the constants,
    %   variables and coefficients in their respective structures. The
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
    %                       - Coefficients for system coefficients that are
    %                           obtained from the parametrization, like
    %                           electric resistance or limiting current
    %                           density. Performing fit using
    %                           fitUI.m function automatically fills the
    %                           values for the coefficients. Alternatively
    %                           they can be manually provided using
    %                           replaceParams method.
    %                       - Dependencies between any two variables. Every
    %                           time any Workspace variables are modified,
    %                           the dependencies are run to refresh the
    %                           dependent variables. Dependencies are
    %                           stored as strings of the commands and
    %                           function calls that are used to calculate
    %                           the dependent variables.
    %
    %   FUNC Methods:
    %       calculate -- Calculates voltage from the UI curve based on
    %                       workspace and user-given values.
    %       copy -- Creates a copy of the object with a new handle.
    %       destructurize -- Destructurizes the structure-based function
    %                           handle to enable fitting with original
    %                           Matlab functions.
    %       removeParams -- Remove parameters from the Workspace structure.
    %       replaceParams -- Replaces the values of existing parameters in
    %                           Workspace structure.
    %       setFuncHandle -- Sets the protected property of funcHandle
    %                           together with the property equation.
    %       setParams -- Sets new parameters to Workspace structure.
    %       viewWorkspace -- Outputs a human-readable report of the
    %                           contents of the Workspace structure.
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
       Workspace; % Structure containing the Coefficients, Variables, Constants and Dependencies between them:
%         Coefficients; % Structure containing the coefficients
%         Variables; % Structure containing the variables
%         Constants; % Structure containing the constants
%         Dependencies; % Structure containing the dependencies as strings
    end
    
    %% Public methods
    methods
        function obj = func(funcHandle,Workspace,varargin)
            % FUNC Constructor method for the object
            %   Inputs:
            %       funcHandle -- The structure-based function handle.
            %       Workspace -- The workspace structure.
            obj.setFuncHandle(funcHandle);
            if func.isWorkspace(Workspace)
                obj.Workspace = Workspace;
            else
                error("Workspace structure of object func should contain exclusively fields 'Coefficients', 'Variables' or 'Constants'.")
            end
            
            obj.refreshWorkspace;
        end
        

        function setFuncHandle(obj,newFuncHandle)
            % SETFUNCHANDLE Sets the funcHandle and equation properties.
            obj.funcHandle = newFuncHandle;
            obj.equation = obj.getEquation;
        end
        
        
        function setParams(obj,SetStruct)
            % SETPARAMS  Sets new parameters in the Workspace structure.
            %  Parameters should be provided as a Workspace-compatible 
            %  structure. The method prioritizes given values over the
            %  existing ones in case of conflict.
            %
            %  See also FUNC.ISWORKSPACE, MERGESTRUCTS
            if func.isWorkspace(SetStruct)
                obj.Workspace = mergeStructs(obj.Workspace,SetStruct);
            else
                error('Parameters to be set have to be given as a single Workspace structure')
            end
            
            obj.refreshWorkspace;
        end
        
        
        function replaceParams(obj,varargin)
            % REPLACEPARAMS  Replaces parameters in the Workspace structure.
            %  Parameters should be provided as name-value pairs or as a
            %  structure (either Workspace compatible or not). Replaces
            %  values only for existing fields with the same name as given.
            %  The method acts recursively inside Workspace.
            %
            % See also ADDVALUESTOSTRUCT
            if isempty(varargin{1})
                return;
            elseif length(varargin) == 1 && isstruct(varargin{1})
                Struct = varargin{1};
                if func.isWorkspace(Struct)
                    fields = fieldnames(Struct);
                    for i = 1:length(fields)
                        replaceParams(obj,Struct.(fields{i}))
                    end
                    return;
                else
                    paramsToReplace = Struct;
                end
            elseif mod(nargin,2)
                for i = 1:2:length(varargin)
                    paramsToReplace.(varargin{i}) = varargin{i+1};
                end
            else
                error('Parameters to be replaced have to be either set as a single structure or as name-value pairs')
            end
            
            if ismember('Dependencies',fieldnames(obj.Workspace))
                TempWorkspace = rmfield(obj.Workspace,'Dependencies'); % Avoid replacing dependencies if named the same as their respective variable
                Dependencies = obj.Workspace.Dependencies;
                TempWorkspace = addValuesToStruct(TempWorkspace,paramsToReplace);
                TempWorkspace.Dependencies = Dependencies;
                obj.Workspace = TempWorkspace;
            else
                obj.Workspace = addValuesToStruct(obj.Workspace,paramsToReplace);
            end
            
            obj.refreshWorkspace;
        end
        
        function removeParams(obj,varargin)
            % REMOVEPARAMS  Removes parameters from the Workspace structure.
            %  Input has to be provided either as a list of
            %  strings/character vectors or as a cell array of
            %  strings/character vectors containing the names of the
            %  parameters to be removed. The method acts recursively inside
            %  the Workspace.
            
            if length(varargin) == 1 && iscell(varargin{1})
                paramsToRemove = varargin{1};
            else
                paramsToRemove = varargin;
            end
            
            for i = 1:length(paramsToRemove)
                if ismember(paramsToRemove{i},fieldnames(obj.Workspace.Constants))
                    obj.Workspace.Constants = rmfield(obj.Workspace.Constants,paramsToRemove{i});
                end
                if ismember(paramsToRemove{i},fieldnames(obj.Workspace.Coefficients))
                    obj.Workspace.Coefficients = rmfield(obj.Workspace.Coefficients,paramsToRemove{i});
                end
                if ismember(paramsToRemove{i},fieldnames(obj.Workspace.Variables))
                    obj.Workspace.Variables = rmfield(obj.Workspace.Variables,paramsToRemove{i});
                end
            end
        end
        
        
        function result = calculate(obj,varargin)
            % CALCULATE Calculates the result from the function handle.
            %   Variables input as name-value pairs are used for the
            %   calculation and any variable required by the function that is
            %   not input to CALCULATE is looked for from the Workspace of the
            %   object.
            
            if ~mod(nargin,2)
                error("Variables have to be given as name-value pairs")
            end
            
            % Check if Workspace was given as an input
            workspaceCall = strcmpi(varargin,'Workspace');
            if any(workspaceCall)
                TempWorkspace = func.createTempWorkspace(varargin{[false workspaceCall(1:end-1)]});
            else
                % Create a temporary workspace that includes the user input
                % variables in addition to the ones found from the func 
                % Workspace.
                TempWorkspace = func.createTempWorkspace(obj.Workspace);
                if nargin>1
                    TempWorkspace = addValuesToStruct(TempWorkspace,...
                        varargin(1:2:end),...
                        varargin(2:2:end));
                end
            end
                        
            % The result is calculated if the TempWorkspace contains all
            % necessary values for the calculation.
            try
                result = obj.funcHandle(TempWorkspace);
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
                coefficientsNamesForFit,...
                problemVariableNames,...
                problemVariables] = destructurize(obj,...
                                                  independentVariableName)
            % DESTRUCTURIZE Destructurizes the structure-based function handle.
            %   Modifies the structure-based function handle to a function
            %   handle with all the constants and scalar variables
            %   calculated in, and the vector variables (problem variables)
            %   and all the fit coefficients used as separate input
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
            %       coefficientNamesForFit -- Lists of the undefined 
            %                                   coefficient names as a cell 
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
            % Coefficients
            coefficientsNamesForFit = {};
            constantCoefficientNames = {};
            constantCoefficients = [];
            if ~isempty(obj.Workspace.Coefficients)
                allCoefficientNames = fieldnames(obj.Workspace.Coefficients);
                for i = 1:length(allCoefficientNames)
                    coefficient = obj.Workspace.Coefficients.(allCoefficientNames{i});
                    if contains(obj.getEquationBody,['Workspace.Coefficients.' allCoefficientNames{i}])
                        % If the coefficient is used in the equation
                        coefficientsNamesForFit = [coefficientsNamesForFit;...
                                                   allCoefficientNames{i}];
                    else % If the coefficient is not used in the equation
                        if ~isempty(coefficient)
                            constantCoefficientNames = [constantCoefficientNames;...
                                                        allCoefficientNames{i}];
                            constantCoefficients = [constantCoefficients;...
                                                    coefficient];
                        else
                            warning("There is a coefficient whose name is not found from the equation but no value has been assigned to it.")
                        end
                    end
                end
            else
                warning("No coefficients defined for the function.")
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

            % Coefficients
            Coefficients = struct();
            for i = 1:length(coefficientsNamesForFit)
                Coefficients.(coefficientsNamesForFit{i}) = sym(coefficientsNamesForFit{i});
            end
            %
            for i = 1:length(constantCoefficientNames)
                Coefficients.(constantCoefficientNames{i}) = constantCoefficients(i);
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
                                                     'Coefficients',Coefficients));
            
            % Modify back to function handle
            varsList = [coefficientsNamesForFit;...
                        problemVariableNames;...
                        independentVariableName];
            destructurizedFuncHandle = matlabFunction(symbolicFunction,...
                                                      'Vars',varsList);
        end
        
        
        
        
        function childFunc = copy(obj,varargin)
            % COPY Creates an independent copy of the object.
            %   Changing properties of the copy doesn't affect the
            %   original.
            %
            % With input parameter 'empty' creates a copy that has no
            %  funcHandle but keeps the Workspace content.
            
            if nargin == 1
                childFunc = func(obj.funcHandle,obj.Workspace);
            elseif strcmp(varargin{1},'empty')
                childFunc = func(@pass,obj.Workspace);
            else
                error('Unknown input')
            end
        end
        
        
        function report = viewWorkspace(obj)
            % VIEWWORKSPACE Outputs a tabular report of the Workspace.
            
            fieldNameConst = fieldnames(obj.Workspace.Constants);
            fieldNameVars = fieldnames(obj.Workspace.Variables);
            fieldNameCoeff = fieldnames(obj.Workspace.Coefficients);
            fields = [strcat(cell(size(fieldNameConst)),'Constants'),fieldNameConst;...
                         strcat(cell(size(fieldNameVars)),'Variables'),fieldNameVars;...
                         strcat(cell(size(fieldNameCoeff)),'Coefficients'),fieldNameCoeff];
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
                    if length(value(:,1)) == 1
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
                else
                    descriptionText = ": no values assigned";
                end
                description{i} = strcat(fields{i,1}(1:end-1),descriptionText);
            end
            report = table(description,valueMean,valueMin,valueMax,standardDeviation);
            report.Properties.RowNames = fieldName;
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
            
            newFunc = func(newFuncHandle,NewWorkspace);
            
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
                'Coefficients',struct()));
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
                             {'Coefficients';...
                             'Variables';...
                             'Constants';...
                             'Dependencies'}));
        end
    end
    
    %% Private dynamic methods
    methods (Access = private)

        function equationStr = getEquation(obj)
            % GETEQUATION Creates human-readable equation string.
            %   Erases all occurences of Workspace from the equation body
            %   leaving a string with human-readable equation structure.            
            equationStr = erase(obj.getEquationBody,...
                                {'Workspace.Coefficients.',...
                                'Workspace.Variables.',...
                                'Workspace.Constants.'});
        end
        
        function equationBody = getEquationBody(obj)
            % GETEQUATIONBODY Outputs only the equation body.
            %  Removes the @(Workspace) from the function handle.            
            equationBody = erase(func2str(obj.funcHandle),'@(Workspace)');
        end
        
        function refreshWorkspace(obj)
            % REFRESHWORKSPACE Recalculates dependent Workspace parameters.
            Workspace = obj.Workspace;
            if ismember('Dependencies',fieldnames(Workspace))
                fn = fieldnames(Workspace.Dependencies);
                for i = 1:numel(fn)
                    try
                        eval(Workspace.Dependencies.(fn{i}))
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:InputParser:ArgumentFailedValidation')||strcmp(ME.identifier,'MATLAB:undefinedVarOrClass')
                            continue;
                        else
                            rethrow(ME)
                        end                           
                    end
                end
                obj.Workspace = Workspace;
            else
                return;
            end
        end
        
    end
    
    %% Private static methods
    
    methods (Access = private, Static)
        function TempWorkspace = createTempWorkspace(Struct)
            % CREATETEMPWORKSPACE Creates a copy of the Workspace structure but without values for standard deviation.            
            fn = fieldnames(Struct);
            if isempty(fn)
                TempWorkspace = struct();
            else
                for i = 1:length(fn)
                    if strcmp(fn{i},'Dependencies')
                        continue;
                    elseif isstruct(Struct.(fn{i}))
                        TempWorkspace.(fn{i}) = func.createTempWorkspace(Struct.(fn{i}));
                    elseif isempty(Struct.(fn{i})) || length(Struct.(fn{i})(1,:)) == 1
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
        end
    end
    
    
end