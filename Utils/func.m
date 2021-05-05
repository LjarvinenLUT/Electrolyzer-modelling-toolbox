classdef func < handle
   
    properties (SetAccess = protected)
       name; % Function handle in string form
       func_handle; % Function handle of the function
       workspace; % Structure containing the function workspace variables
    end
    
    methods
        function obj = func(func_handle,workspace)
            obj.func_handle = func_handle;
            obj.func_name = fun2str(obj.func_handle);
            obj.workspace = workspace;
        end
        
        function result = calculate(obj)
            result = obj.func_handle(obj.workspace);
        end
        
    end
end