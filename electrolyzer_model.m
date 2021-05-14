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
       potentialFunc; % Func object for combined overpotential function
%        coeffs; % Structure for model coefficients.
%        vars;   % System measured variables
%        potentials; % Structure array of all potentials
    end
    
    
    methods
       
        % Constructor function
        function obj = electrolyzer_model(varargin)
            defaultType = "pem";
            defaultElectrolyte = "KOH";
            
            parser = inputParser;
            addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
            addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
            
            parse(parser,varargin{:});
            
            parser.Results.type;
            set_type(obj, parser.Results.type);
            set_electrolyte(obj, string(parser.Results.electrolyte));
            
            obj.potentialFunc = "Potential function not defined";
%             obj.potentials = struct('name',{},'func',{},'coeffs',{});
%             obj.coeffs = struct([]);
%             if strcmp(obj.type,"pem")
%                 obj.vars = struct('T',[],'Current',[],'Voltage',[],'p1',[],'p2',[]);
%             elseif strcmp(obj.type,"alkaline")
%                 obj.vars = struct('T',[],'Current',[],'Voltage',[],'p',[],'m',[]);
%             else
%                 obj.vars = struct('T',[],'Current',[],'Voltage',[]);
%             end
        end
        
        % Function for setting electrolyzer type
        function set_type(obj,type)
            obj.type = string(lower(type));
        end
        
        % Function for setting electrolyte for alkali electrolyzers
        function set_electrolyte(obj,electrolyte)
            if strcmp(obj.type,"alkaline")
                obj.electrolyte = electrolyte;
            else
                obj.electrolyte = "none";
            end
        end
        
        % Function for ssetting measured variables
        function set_vars(obj,vars)
            obj.vars = vars;
        end
        
        function add_potential(obj, addedPotentialFunc)
%             obj.potentials = [obj.potentials,added_potential];
%             obj.coeffs = mergeStructs(added_potential.coeffs,obj.coeffs);
            if ~isa(addedPotentialFunc,'func')
                error("Potential to be added has to be in a form of a func object")
            elseif ~isa(obj.potentialFunc,'func')
                obj.potentialFunc = addedPotentialFunc;
            else
                obj.potentialFunc = addFuncs(obj.potentialFunc,addedPotentialFunc);
            end
        end

        
        % Function for fitting UI curve
        function fit_UI(obj,U,I,varargin)
            
            defaultMethod = 'PS';
            
            parser = inputParser;
            addRequired(parser,'obj',@(x) isa(x,'electrolyzer_model'))
            addRequired(parser,'U',@(x) isnumeric(x))
            addRequired(parser,'I',@(x) isnumeric(x))
            addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x)) % Fitting method to be used
            
            parse(parser,obj,U,I,varargin{:});
            
            method = upper(string(parser.Results.method));
            
            obj.fit_parameters = fit_UI(obj.potentialFunc,U,I,'method',method);
            
        end
        
        % Function for plotting UI curve
        function show_UI(obj)
            
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
        function clearOverpotentials(obj)
           obj.potentialFunc = "Potential function not defined"; 
        end
        
    end
    
end