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
%           TODO: Check that no field in the Workspace structures is empty
            result = obj.funcHandle(obj.Workspace);
        end
        
        function setFuncHandle(obj,newFuncHandle)
            obj.funcHandle = newFuncHandle;
            obj.equation = func2str(obj.funcHandle);
        end
        
        function equationStr = getEquation(obj)
            equationStr = erase(obj.getEquationBody,{'Workspace.Coefficients.','Workspace.Variables.','Workspace.Constants.'});
        end
        
        function equationBody = getEquationBody(obj)
            equationBody = erase(obj.equation,'@(Workspace)');
        end
        
    end
end