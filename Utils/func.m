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
    %   inputs.
    %
    %   FUNC Properties:
    %       equation -- Function handle as a easily human-readable string.
    %       funcHandle -- The structure based function handle that uses
    %                       only Workspace as an input.
    %       Workspace -- The workspace structure containing substructures:
    %                       - Constants for system related constants,
    %                       - Variables for measured variable values,
    %                       - Coefficients for system coefficients.
    %
    %   FUNC Methods:
    %       calculate -- Calculates voltage from the UI curve based on
    %                       workspace and user-given values.
    %       copy -- Creates a copy of the object with a new handle.
    %       destructurize -- Destructurizes the structure-based function
    %                           handle to enable fitting with original
    %                           Matlab functions.
    %       getEquationBody -- Returns a string containing only the body of
    %                           the function handle by removing
    %                           @(Workspace).
    %       replaceParams -- Replaces the values of existing parameters in
    %                           Workspace structure
    %       setFuncHandle -- Sets the protected property of funcHandle
    %                           together with the property equation.
    %       setParams -- Sets new parameters to Workspace structure
    %       viewWorkspace -- Outputs a human-readable raport of the
    %                           contents of the Workspace structure.
    %   FUNC static methods:
    %       add -- Use addition to combine two func objects into a new one.
    %
    %   See also FUNCTION_HANDLE, ADDFUNCS, ISCOMPLETESTRUCT
   
    properties (SetAccess = protected)
       equation; % Function handle in string form
       funcHandle; % Function handle of the function that uses only Workspace structure as an input
       Workspace; % Structure containing the Coefficients, Variables and Constants:
%         Coefficients; % Structure containing the coefficients
%         Variables; % Structure containing the variables
%         Constants; % Structure containing the constants
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
        end
        

        function setFuncHandle(obj,newFuncHandle)
            % SETFUNCHANDLE Sets the property of funcHandle together with the property of equation.
            obj.funcHandle = newFuncHandle;
            obj.equation = obj.getEquation;
        end
        
        
        function equationBody = getEquationBody(obj)
            % GETEQUATIONBODY Outputs the equation body from the function handle without @(Workspace).
            equationBody = erase(func2str(obj.funcHandle),'@(Workspace)');
        end
        
        function replaceParams(obj,varargin)
            % REPLACEPARAMS  A method for replacing parameters in the Workspace structure.
            %  Parameters should be provided as a name-value pair or
            %  directly as a structure.
            if isempty(varargin{1})
                return;
            elseif length(varargin) == 1 && isstruct(varargin{1})
                paramsToReplace = varargin{1};
            elseif mod(nargin,2)
                for i = 1:2:length(varargin)
                    paramsToReplace.(varargin{i}) = varargin{i+1};
                end
            else
                error('Parameters to be replaced have to be either set as a single structure or as name-value pairs')
            end
            
            obj.Workspace = addValuesToStruct(obj.Workspace,paramsToReplace);
        end
        
        function removeParams(obj,varargin)
            % REMOVEPARAMS  A method for removing parameters from the Workspace structure.
            %  Parameter names should be provided as a list or a cell array
            %  of strings.
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
        
        function setParams(obj,SetStruct)
            % SETPARAMS  A method for setting parameters in the Workspace structure.
            %  Parameters should be provided as a structure.
            if func.isWorkspace(SetStruct)
                obj.Workspace = mergeStructs(obj.Workspace,SetStruct);
            else
                error('Parameters to be set have to be given as a single Workspace structure')
            end
        end
        
        
        function result = calculate(obj,varargin)
            % CALCULATE Calculates the voltage from the UI curve.
            %   Variables input as name-value pairs are used for the
            %   calculation and any variable required by the function that is
            %   not input to CALCULATE is looked for from the Workspace of the
            %   object.
            
            if ~mod(nargin,2)
                error("Variables have to be given as name-value pairs")
            end
            
            % Create a temporary workspace that includes the user input 
            % variables in addition to the ones found from.
            TempWorkspace = func.createTempWorkspace(obj.Workspace);
            if nargin>1
                TempWorkspace = addValuesToStruct(TempWorkspace,...
                                                  varargin(1:2:end),...
                                                  varargin(2:2:end));
            end
            
            % The result is calculated if the TempWorkspace contains all
            % necessary values in all the fields.
            if isCompleteStruct(TempWorkspace)
                result = obj.funcHandle(TempWorkspace);
            else
                error('Result cannot be calculated because one or more workspace entries are missing. Either add the missing entries to the workspace of the func object or input them to calculate method as name-value pairs.')
            end
        end
        
        
        function [destructurizedFuncHandle,...
                coefficientsNamesForFit,...
                problemVariableNames,...
                problemVariables] = destructurize(obj,...
                                                  independentVariableName)
            % DESTRUCTURIZE Destructurizes the structure-based function handle.
            %   Modifies the function handle to use separate input
            %   parameters instead of a  singleWorkspace structure.
            %   Destructurizing enable fitting with original Matlab 
            %   functions like fit and particleswarm. 
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
                if strcmp(allVariableNames{i},independentVariableName)
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
            % COPY Creates a copy of the object with a new handle. With
            %   parameter 'empty' create a copy that has no funcHandle but
            %   keeps the Workspace content.
            if nargin == 1
                childFunc = func(obj.funcHandle,obj.Workspace);
            elseif strcmp(varargin{1},'empty')
                childFunc = func(@pass,obj.Workspace);
            else
                error('Unknown input')
            end
        end
        
        
        function report = viewWorkspace(obj)
            % VIEWWORKSPACE Outputs a tabular report of the contents of the Workspace structure.
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
            % Constants
            for i = 1:length(fieldName)
                value = obj.Workspace.(fields{i,1}).(fields{i,2});
                if ~isempty(value)
                    switch length(value(:,1))
                        case 1
                            valueMean(i) = value(:,1);
                            descriptionText = ": scalar";
                            if length(value(1,:)) == 2
                                standardDeviation(i) = value(1,2);
                                descriptionText = strcat(descriptionText," with confidence bounds");
                            else
                                descriptionText = strcat(descriptionText," without confidence bounds");
                            end
                        otherwise
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
    
    %% Private dynamic methods
    methods (Access = private)

        function equationStr = getEquation(obj)
            % GETEQUATION Erases all occurences of Workspace from the
            %   equation body leaving a string with human-readable
            %   equation structure.
            equationStr = erase(obj.getEquationBody,...
                                {'Workspace.Coefficients.',...
                                'Workspace.Variables.',...
                                'Workspace.Constants.'});
        end
    end
    
    %% Static methods
    
    methods(Static, Access = public)
        function newFunc = add(func1,func2)
            % ADD Add two func objects together combining their workspaces and
            %   using addition for combining their function handles.
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
            
        end
        
        function emptyFunc = createEmpty
            % CREATEEMPTY Create an empty func object with a function
            % handle that does nothing.
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
            % ISWORKSPACE Evaluates if a given structure fulfills
            %   requirements for a Workspace.
            b = all(ismember(fieldnames(Struct),...
                             {'Coefficients';...
                             'Variables';...
                             'Constants'}));
        end
    end
    
    methods (Access = private, Static)
        function TempWorkspace = createTempWorkspace(Struct)
            % CREATETEMPWORKSPACE Creates a copy of the Workspace structure
            %   but without values for standard deviation
            fn = fieldnames(Struct);
            if isempty(fn)
                TempWorkspace = struct();
            else
                for i = 1:length(fn)
                    if isstruct(Struct.(fn{i}))
                        TempWorkspace.(fn{i}) = func.createTempWorkspace(Struct.(fn{i}));
                    elseif isempty(Struct.(fn{i})) || length(Struct.(fn{i})(1,:)) == 1
                        TempWorkspace.(fn{i}) = Struct.(fn{i});
                    else
                        TempWorkspace.(fn{i}) = Struct.(fn{i})(:,1);
                    end
                end
            end
        end
    end
    
    
end