% Matlab programming style guideline could be used when code is ready to
% format everything to standard format
% Guide can be found in Zotero folder "Electrolyzer modeling/matlab style guide"

%% Class for electrolyzer modelling
% (This should probably be just a superclass and subclasses for each type should be separate)

% Make the class a handle class so that variables can be updated
classdef electrolyzer_model < handle
    
    properties (SetAccess = protected)
       type; % PEM or alkaline
       electrolyte; % Electrolyte for alkali electrolyzers (should be in its own subclass)
       overpotential_function; % Function handle for combined overpotential function
       fit_parameters; % Table for fitted parameters. 
                    % Using table for parameters enables use of
                    % models with varying number of fit parameters and
                    % their different namings
       T;   % System temperature
       
       overpotentials; % Array of all overpotential function handles
    end
    
    
    methods
       
        % Constructor function
        function obj = electrolyzer_model(varargin)
            defaultElectrolyte = "KOH";
            
            parser = inputParser;
            addParameter(parser,'type',@(x) ischar(x)||isstring(x));
            addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
            
            parse(parser,varargin{:});
            
            parser.Results.type
            set_type(obj, parser.Results.type);
            set_electrolyte(obj, string(parser.Results.electrolyte));
        end
        
        % Function for setting electrolyzer type
        function set_type(obj,type)
            obj.type = string(lower(type));
        end
        
        % Function for setting electrolyte for alkali electrolyzers
        function set_electrolyte(obj,electrolyte)
            if ~strcmp(obj.type,"alkaline")
                error("Setting electrolyte possible only with alkali electrolyzers")
            end
            obj.electrolyte = electrolyte;
        end
        
        function object = add_overpotential(obj, overpotential)
            obj.overpotentials{end+1} = overpotential;
            object = obj;
        end
        
        % Function for combining overpotential function handles
        function object = combine_overpotentials(obj)
            % Call utility function to combine overpotential function
            % handles to on function
            obj.overpotential_function = combineFuncHandles(obj.overpotentials);
            object = obj;
        end
        
        % Function for fitting UI curve
        function object = fit_UI(obj,U,I,varargin)
            
            defaultMethod = 'PS';
            
            parser = inputParser;
            addRequired(parser,'obj',@(x) isa(x,'electrolyzer_model'))
            addRequired(parser,'U',@(x) isnumeric(x))
            addRequired(parser,'I',@(x) isnumeric(x))
            addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x)) % Fitting method to be used
            
            parse(parser,obj,U,I,varargin{:});
            
            method = upper(string(parser.Results.method));
            
            obj.fit_parameters = fit_UI(obj.overpotential_function,U,I,'method',method);
            
            object = obj;
        end
        
        % Function for plotting UI curve
        function object = show_UI(obj)
            display("Print figures")
            
            object = obj;
        end
        
        % Get all unique arguments from the given overpotential function
        % handles
        function uniqueArguments = getOverpotentialArguments(obj)
            addpath("Utils");
            argumentMap = containers.Map();
            for i = 1:length(obj.overpotentials)
                arguments = getFunctionArguments(obj.overpotentials{i});
                
                % Loop through all the arguments
                for j = 1:length(arguments)
                    % By using map we dont have to check if the argument
                    % has already been added
                    argumentMap(arguments{j}) = arguments{j};
                end
            end
            uniqueArguments = values(argumentMap);
        end
        
        % Clears the overpotential cell array
        function clear_overpotentials(obj)
           obj.overpotentials = {}; 
        end
        
    end
    
end