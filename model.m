% Electrolyzer model superclass

classdef model < handle
    
    properties
        Variables; % Table of variables (temperature, current, voltage, pressure) that have been measured from the system
        Parameters; % Table of model parameters
    end
    
    properties (SetAccess = private)
        U; % Modelled voltage vector
    end
    
    
    methods
       
        % Model initialiser
        function obj = model(V)
        
            obj.Variables = V;
        
        end
        
        % Modelled voltage calculation
        function U = 
        
        
        
    end
    
end