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
        
        function result = calculate(obj)
            if isCompleteStructure(obj.Workspace)
                result = obj.funcHandle(obj.Workspace);
            else
                warning('Result cannot be calculated because one or more workspace entries are missing.')
                result = nan;
            end
        end
        
        function setFuncHandle(obj,newFuncHandle)
            obj.funcHandle = newFuncHandle;
            obj.equation = obj.getEquation;
        end
        
        function equationStr = getEquation(obj)
            equationStr = erase(obj.getEquationBody,{'Workspace.Coefficients.','Workspace.Variables.','Workspace.Constants.'});
        end
        
        function equationBody = getEquationBody(obj)
            equationBody = erase(func2str(obj.funcHandle),'@(Workspace)');
        end
        
        function [destructurizedFuncHandle,coefficientsNamesForFit,independentVariable,dependentVariable,problemVariableNames,problemVariables] = destructurizeFuncHandle(obj,independentVariableName,dependentVariableName)
            % Constants
            if ~isCompleteStructure(obj.Workspace.Constants)
                error('Missing constant values. Not able to destructurize function handle for fitting.')
            end
            
            % Coefficients
            coefficientsNamesForFit = {};
            allCoefficientNames = fieldnames(obj.Workspace.Coefficients);
            for i = 1:length(allCoefficientNames)
                if isempty(obj.Workspace.Coefficients.(allCoefficientNames{i}))
                    coefficientsNamesForFit = [coefficientsNamesForFit;allCoefficientNames{i}];
                end
            end
            
            % Variables
            allVariableNames = fieldnames(obj.Workspace.Variables);
            if ~ismember(dependentVariableName,allVariableNames)
                error('Dependent variable is not included in the workspace. Not able to destructurize function handle for fitting.')
            elseif isempty(obj.Workspace.Variables.(dependentVariableName))
                error("Value of the dependent variable, '" + dependentVariableName + "', is not defined. Not able to destructurize function handle for fitting.")
            else
                dependentVariable = obj.Workspace.Variables.(dependentVariableName);
            end
            if ~ismember(independentVariableName,allVariableNames)
                error('Inependent variable is not included in the workspace. Not able to destructurize function handle for fitting.')
            elseif isempty(obj.Workspace.Variables.(independentVariableName))
                error("Value of the independent variable, '" + independentVariableName + "', is not defined. Not able to destructurize function handle for fitting.")
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
                if any(strcmp(allVariableNames{i},{dependentVariableName,independentVariableName}))
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
            
            % Create destructurized function handle
            
        end
        
    end
end