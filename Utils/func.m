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
    %       setFuncHandle -- Sets the protected property of funcHandle
    %                           together with the property equation.
    %       viewWorkspace -- Outputs a human-readable raport of the
    %                           contents of the Workspace structure.
    %
    %   See also ADDFUNCS
   
    properties (SetAccess = protected)
       equation; % Function handle in string form
       funcHandle; % Function handle of the function that uses only Workspace structure as an input
    end
    
    properties
       Workspace; % Structure containing the Coefficients, Variables and Constants:
%         Coefficients; % Structure containing the coefficients
%         Variables; % Structure containing the variables
%         Constants; % Structure containing the constants
    end
    
    %% Public methods
    methods
        function obj = func(funcHandle,Workspace)
            % FUNC Constructor method for the object
            %   Inputs:
            %       funcHandle -- The structure-based function handle.
            %       Workspace -- The workspace structure.
            obj.setFuncHandle(funcHandle);
            if any(~ismember(fieldnames(Workspace),{'Coefficients';'Variables';'Constants'}))
                error("Workspace structure of object func should contain exclusively fields 'Coefficients', 'Variables' or 'Constants'.")
            else
                obj.Workspace = Workspace;
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
            TempWorkspace = obj.Workspace;
            if nargin>1
                TempWorkspace = addValuesToStruct(TempWorkspace,varargin(1:2:end),varargin(2:2:end));
            end
            
            % The result is calculated if the TempWorkspace contains all
            % necessary values in all the fields.
            if isCompleteStructure(TempWorkspace)
                result = obj.funcHandle(TempWorkspace);
            else
                error('Result cannot be calculated because one or more workspace entries are missing. Either add the missing entries to the workspace of the func object or input them to calculate method as name-value pairs.')
            end
        end
        
        
        function [destructurizedFuncHandle,coefficientsNamesForFit,problemVariableNames,problemVariables] = destructurize(obj,independentVariableName)
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
            if ~isCompleteStructure(obj.Workspace.Constants)
                error('Missing constant values. Not able to destructurize function handle for fitting.')
            end
            %
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

            % Create destructurized function handle

            % Constants
            Constants = obj.Workspace.Constants;

            % Coefficients
            Coefficients = struct();
            for i = 1:length(coefficientsNamesForFit)
                Coefficients.(coefficientsNamesForFit{i}) = sym(coefficientsNamesForFit{i});
            end
            %
            for i = 1:length(predefinedCoefficientNames)
                Coefficients.(predefinedCoefficientNames{i}) = predefinedCoefficients(i);
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
            symbolicFunction = obj.funcHandle(struct('Constants',Constants,'Variables',Variables,'Coefficients',Coefficients));
            
            % Modify back to function handle
            varsList = [coefficientsNamesForFit;problemVariableNames;independentVariableName];
            destructurizedFuncHandle = matlabFunction(symbolicFunction,'Vars',varsList);
        end
        

        function childFunc = copy(obj)
            % COPY Creates a copy of the object with a new handle.
            childFunc = func(obj.funcHandle,obj.Workspace);
        end
        
        
        function output = viewWorkspace(obj)
            % VIEWWORKSPACE Outputs a human-readable raport of the
            %   contents of the Workspace structure.
            % TODO
        end
        
    end
    
    %% Private methods
    methods (Access = private)

        function equationStr = getEquation(obj)
            % GETEQUATION Erases all occurences of Workspace from the
            %   equation body leaving a string with human-readable
            %   equation structure.
            equationStr = erase(obj.getEquationBody,{'Workspace.Coefficients.','Workspace.Variables.','Workspace.Constants.'});
        end
    end
end