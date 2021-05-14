classdef func < handle
   
    properties (SetAccess = protected)
       equation; % Function handle in string form
       funcHandle; % Function handle of the function that uses only Workspace structure as an input
    end
    
    properties
       Workspace; % Structure containing:
%         Coefficients; % Structure containing the coefficients
%         Variables; % Structure containing the variables
%         Constants; % Structure containing the function workspace variables
    end
    
    methods
        function obj = func(funcHandle,Workspace)
            obj.setFuncHandle(funcHandle);
            if any(~ismember(fieldnames(Workspace),{'Coefficients';'Variables';'Constants'}))
                error("Workspace structure of object func should contain exclusively fields 'Coefficients', 'Variables' or 'Constants'.")
            else
                obj.Workspace = Workspace;
            end
        end
        
       
        
        %%
        function setFuncHandle(obj,newFuncHandle)
            obj.funcHandle = newFuncHandle;
            obj.equation = obj.getEquation;
        end
        

        %%
        function equationStr = getEquation(obj)
            equationStr = erase(obj.getEquationBody,{'Workspace.Coefficients.','Workspace.Variables.','Workspace.Constants.'});
        end
        
        %%
        function equationBody = getEquationBody(obj)
            equationBody = erase(func2str(obj.funcHandle),'@(Workspace)');
        end
        
         %%
        function result = calculate(obj,varargin)
            
            if ~mod(nargin,2)
                error("Variables have to be given as name-value pairs")
            end
            
            TempWorkspace = obj.Workspace;
            if nargin>1
                TempWorkspace = addValuesToStruct(TempWorkspace,varargin(1:2:end),varargin(2:2:end));
            end
                
            if isCompleteStructure(TempWorkspace)
                result = obj.funcHandle(TempWorkspace);
            else
                error('Result cannot be calculated because one or more workspace entries are missing. Either add the missing entries to the workspace of the func object or input them to calculate method as name-value pairs.')
            end
        end
        
        
        %%
        function [destructurizedFuncHandle,coefficientsNamesForFit,problemVariableNames,problemVariables] = destructurize(obj,independentVariableName)
            % Constants
            if ~isCompleteStructure(obj.Workspace.Constants)
                error('Missing constant values. Not able to destructurize function handle for fitting.')
            end
            
            % Coefficients
            coefficientsNamesForFit = {};
            predefinedCoefficientNames = {};
            predefinedCoefficients = [];
            predefCoeffsExist = false;
            if ~isempty(obj.Workspace.Coefficients)
                allCoefficientNames = fieldnames(obj.Workspace.Coefficients);
                for i = 1:length(allCoefficientNames)
                    coefficient = obj.Workspace.Coefficients.(allCoefficientNames{i});
                    if isempty(coefficient)
                        coefficientsNamesForFit = [coefficientsNamesForFit;allCoefficientNames{i}];
                    else
                        predefCoeffsExist = true;
                        predefinedCoefficientNames = [predefinedCoefficientNames;allCoefficientNames{i}];
                        predefinedCoefficients = [predefinedCoefficients;coefficient];
                    end
                end
            end
            
            if predefCoeffsExist
                warning("There were one or more predefined coefficients in the workspace. The predefined value for these coefficients have been used.")
            end
            
            
            % Variables
            allVariableNames = fieldnames(obj.Workspace.Variables);
            if ~ismember(independentVariableName,allVariableNames)
                error('Inependent variable is not included in the workspace. Not able to destructurize function handle for fitting.')
            else
                independentVariable = obj.Workspace.Variables.(independentVariableName);
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
                    nonProblemVariableNames = [nonProblemVariableNames;variableName];
                    nonProblemVariables = [nonProblemVariables;variable];
                else
                    problemVariableNames = [problemVariableNames;variableName];
                    problemVariables = [problemVariables;variable];
                end                    
            end
            
            
            %% Create destructurized function handle
            
            % Constants
            Constants = obj.Workspace.Constants;
            
            % Coefficients
            Coefficients = struct();
            for i = 1:length(coefficientsNamesForFit)
                Coefficients.(coefficientsNamesForFit{i}) = sym(coefficientsNamesForFit{i});
            end
            
            for i = 1:length(predefinedCoefficientNames)
                Coefficients.(predefinedCoefficientNames{i}) = predefinedCoefficients(i);
            end
            
            % Variables
            for i = 1:length(nonProblemVariableNames)
                Variables.(nonProblemVariableNames{i}) = nonProblemVariables(i);
            end
            
            for i = 1:length(problemVariableNames)
                Variables.(problemVariableNames{i}) = sym(problemVariableNames{i});
            end
            
            Variables.(independentVariableName) = sym(independentVariableName);
            
            % Build new function handle in symbolic form
            symbolicFunction = obj.funcHandle(struct('Constants',Constants,'Variables',Variables,'Coefficients',Coefficients));
            
            % Modify back to function handle
            varsList = [coefficientsNamesForFit;problemVariableNames;independentVariableName];
            destructurizedFuncHandle = matlabFunction(symbolicFunction,'Vars',varsList);
            
        end
        
        %%
        function childFunc = copy(obj)
            childFunc = func(obj.funcHandle,obj.Workspace);
        end
        
        %%
        function output = viewWorkspace(obj)
            % TODO
        end
        
    end
end