%% Class for electrolyzer modelling
% (This should probably be just a superclass and subclasses for each type should be separate)

% Make the class a handle class so that variables can be updated
classdef electrolyzer_model < handle
    
    properties (SetAccess = protected)
       type; % PEM or alkaline
       electrolyte; % Electrolyte for alkali electrolyzers (should be in its own subclass)
       Overpotential_function; % Function handle for overpotentials
       parameters; % Structure array for fitted parameters. 
                    % Using structure array for parameters enables use of
                    % models with varying number of fit parameters and
                    % their different namings
    end
    
    
    methods
       
        % Creator function
        function obj = electrolyzer_model(varargin)

            parser = inputParser;
            addParameter(parser,'type',@(x) ischar(x)||isstring(x));
            addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
            
            parse(parser,varargin{:});
            
            
            set_type(parser.Results.type);
            set_electrolyte(string(parser.Results.electrolyte));
        end
        
        % Function for setting electrolyzer type
        function set_type(obj,type)
            obj.type = string(lower(type));
        end
        
        % Function for setting electrolyte for alkali electrolyzers
        function set_electrolyte(obj,electrolyte)
            if ~strcmp(obj.type,"alkali")
                error("Setting electrolyte possible only with alkali electrolyzers")
            end
            obj.electrolyte = electrolyte;
        end
        
        % Function for setting overpotential function handles
        function set_overpotentials(obj)
            
        end
        
        % Function for fitting UI curve
        function fit_UI(obj,U,I)
            
        end
        
        % Function for plotting UI curve
        function plot_UI(obj)
            
        end
        
    end
    
end