classdef electrolyzerModel < handle
    % ELECTROLYZERMODEL Contains all the necessary information about the electrolyzer model used and enables its simple use.
    %
    %   An object of ELECTROLYZERMODEL can be used to store measured
    %   variables and the functions determining the UI curve.
    %   ELECTROLYZERMODEL also enables fitting of UI curve parameters and
    %   calculating the voltage from input current based on the UI curve.
    %
    %   ELECTROLYZERMODEL Properties:
    %       electrolyte -- Chemical composition of the used electrolyte.
    %                       Significant to alkali electrolysis only.
    %       potentialFunc -- A func object containing the UI curve
    %                           determining equation, values of its
    %                           variables, constants and coefficients, and
    %                           functionality for UI fitting and
    %                           calculation.
    %       type -- Type of the electrolyzer, PEM or alkaline
    %       Variables -- A structure containing measured variables for the
    %                       model
    %
    %   ELECTROLYZERMODEL Methods:
    %       addAllPotentials -- Adds all the potentials listed as input
    %                           to the potentialFunc object with their
    %                           default settings.
    %       addPotential -- Adds an expression of a given potential to
    %                       the potentialFunc object.
    %       calculate -- Calculates voltage from the UI curve for given 
	%                       currents.
    %       clearOverpotentials -- Clears the potentialFunc property.
    %       copy -- Creates a copy of the object whose properties are no
    %               longer linked to the parent.
    %       fitUI -- Performs a fit for the UI curve extracting values for
    %                   its coefficients.
    %       removeVariables -- Remove variables from the Variables
    %                           structure.
    %       showCoefficients -- Lists fit coefficient values.
    %       showUI -- Plots the UI curve for the model.
    %       setVariables -- Set variables in the Variables structure
    %
    %   See also FUNC
    
    %% Properties
    properties (SetAccess = protected)
       type; % PEM or alkaline
       electrolyte; % Electrolyte
       potentialFunc; % A func object for combined overpotential function
       Variables; % A structure that contains the measured variables for the model
    end
    
    %% Public methods
    methods
        
        function obj = electrolyzerModel(varargin)
%           ELECTROLYZERMODEL  Constructor function for the object.
%               Optional inputs:
%                   type -- electrolyzer type.
%                               Default: PEM
%                   electrolyte -- Chemical composition of the electrolyte.
%                                   Default: KOH
            defaultType = "pem";
            defaultElectrolyte = "KOH";
            
            parser = inputParser;
            addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
            addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
            
            parse(parser,varargin{:});
            
            parser.Results.type;
            setType(obj, parser.Results.type);
            setElectrolyte(obj, string(parser.Results.electrolyte));
            
            obj.potentialFunc = "Potential function not defined";
        end
        
        function setVariables(obj,varargin)
%            SETVARIABLES  A method for setting variables in the Variables structure. 
%               Variables should be provided as a name-value pair or 
%               directly as a structure.
            if isempty(varargin{1})
                return;
            elseif length(varargin) == 1 && isstruct(varargin{1})
                obj.Variables = varargin{1};
            elseif mod(nargin,2)
                for i = 1:2:length(varargin)
                    obj.Variables.(varargin{i}) = varargin{i+1};
                end
            else
                error('Variables have to be either set as a single structure or as name-value pairs')
            end
        end
        
        function removeVariables(obj,varargin)
%           REMOVEVARIABLES  A method for removing variables from the Variables structure. 
%               Variables to be removed have to be provided as strings
            if length(varargin) == 1 && iscell(varargin{1})
                variablesToRemove = varargin{1};
            else
                variablesToRemove = varargin;
            end
            
            obj.Variables = rmfield(obj.Variables,variablesToRemove);
        end
        
        
        function showCoefficients(obj)
%            SHOWCOEFFICIENTS  A method for showing fit UI curve coefficients. 
%               (TODO) coefficients found from the workspace of potentialFunc
        end
        
        
        function addPotential(obj, argin)
%           ADDPOTENTIAL  Adds the given potential term to the total potential func object. 
%               Input of a string uses getPotential function to get the 
%               default func object. Alternatively the user can input a 
%               func object directly.
%
%           Examples:
%
%           obj.ADDPOTENTIAL('nernst') adds the Nernst potential term to
%               the potentialFunc parameter using variables and 
%               electrolyzer type defined for the electrolyzerModel object.
%
%           obj.ADDPOTENTIAL(func) adds the potential term defined by the
%               input func object to potentialFunc parameter.

            if isstring(argin) || ischar(argin)
                addedPotentialFunc = obj.getPotential(argin);
            elseif isa(argin,'func')
                addedPotentialFunc = argin;
            else
                error("Potential to be added has to be a func object, or you have to specify with a string which potential term you want to add")
            end
            
            if ~isa(obj.potentialFunc,'func')
                obj.potentialFunc = addedPotentialFunc;
            else
                obj.potentialFunc = addFuncs(obj.potentialFunc,addedPotentialFunc);
            end
        end
        
        function addAllPotentials(obj, varargin)
%            ADDALLPOTENTIALS Adds all the potentials that user has requested.
%               Names to be added are provided in an input cell array or 
%               as separate input values.
            if length(varargin) == 1 && iscell(varargin{1})
                potentialNames = varargin{1};
            else
                potentialNames = varargin;
            end
            for i = 1:length(potentialNames)
                obj.addPotential(potentialNames{i})
            end
        end

        
        function [fitCoefficients,gof] = fitUI(obj,varargin)
%           FITUI A method for extracting the fit coefficients for the electrolyzerModel.
%               Calls external function fitUI.
%               
%               obj.FITUI() performs UI curve fit using voltage and current
%                   data contained in the Variables structure of the
%                   electrolyzerModel object. Default method of particle
%                   swarm optimisation of the sum of square residuals is 
%                   used for the fitting.
%
%               obj.FITUI(U,I) performs UI curve fit using voltage and
%                   current data provided as input parameters.
%
%               obj.FITUI(_,'method',m) performs UI durve fit using the
%                   method defined as input m. 
%                   Options:
%                       PS -- Particle swarm optimisation
%                       NLLSE -- Non linear least squares error regression
            
            defaultMethod = 'PS';
            defaultWeights = 'default';
            defaultU = nan;
            defaultI = nan;
            
            parser = inputParser;
            addRequired(parser,'obj',@(x) isa(x,'electrolyzerModel'))
            addOptional(parser,'U',defaultU,@(x) isnumeric(x))
            addOptional(parser,'I',defaultI,@(x) isnumeric(x))
            addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x)) % Fitting method to be used
            addParameter(parser,'weights',defaultWeights,@(x) ischar(x)||isstring(x))

            parse(parser,obj,varargin{:});
            
            U = parser.Results.U;
            I = parser.Results.I;
            method = upper(string(parser.Results.method));
            weightsMethod = lower(string(parser.Results.weights));
            
            if isnan(I) % Current not provided as an input
                INameOptions = {'I','current','j','i'};
                memberNameIndexI = ismember(INameOptions,fieldnames(obj.Variables));
                nMemberNameIndexI = sum(memberNameIndexI);
                switch nMemberNameIndexI
                    case 0
                        error("Current not defined. Defined it either in the Variables structure of the electrolyzerModel object. Alternatively, if voltage has been provided as the second parameter for fitUI, current can be provided as the third parameter.")
                    case 1
                        I = obj.Variables.(INameOptions{memberNameIndexI});
                    otherwise
                        error("Current overdefined in the Variables structure of the electrolyzerModel object.")
                end
                if isnan(U) % Voltage not provided as an input
                    UNameOptions = {'U','voltage','V','u'};
                    memberNameIndexU = ismember(UNameOptions,fieldnames(obj.Variables));
                    nMemberNameIndexU = sum(memberNameIndexU);
                    switch nMemberNameIndexU
                        case 0
                            error("Voltage not defined. Defined it either in the Variables structure of the electrolyzerModel object. Alternatively, if voltage has been provided as the second parameter for fitUI, current can be provided as the third parameter.")
                        case 1
                            U = obj.Variables.(UNameOptions{memberNameIndexU});
                        otherwise
                            error("Voltage overdefined in the Variables structure of the electrolyzerModel object.")
                    end
                end
            end
            
            [fitCoefficients,gof] = fitUI(obj.potentialFunc,U,I,'method',method,'weights',weightsMethod);
            
            
        end
        
        function result = calculate(obj,varargin)
            % CALCULATE Calculates voltage from the UI curve.
            %   Calls the calculate method of the func object
            %   potentialFunc.
            %   Variables input as name-value pairs are used for the
            %   calculation and any variable required by the function that is
            %   not input to CALCULATE is looked for from the Workspace of the
            %   potentialFunc.
            result = obj.potentialFunc.calculate(varargin{:});
        end
        
        
        function showUI(obj)
            % SHOWUI (TODO) Creates a figure and plots the UI curve on it
            
        end
        
        
        function clearOverpotentials(obj)
            % CLEAROVERPOTENTIALS Clears the potential function
            obj.potentialFunc = "Potential function not defined";
        end
        
        function childObj = copy(obj)
            % COPY  Creates a full copy of the object with its own handle.
            %   Properties of the child object are no longer related to
            %   those of the parent.
            childObj = electrolyzerModel('type',obj.type,'electrolyte',obj.electrolyte);
            childObj.addPotential(obj.potentialFunc.copy);
            childObj.setVariables(obj.Variables);
        end
    end
    
    %% Private methods
    methods (Access = private)
        
        function setType(obj,type)
%           SETTYPE  Function for setting electrolyzer type

            if strcmpi(type,'alkali')
                type = 'alkaline';
            elseif ~any(strcmpi(type,{'pem','alkaline'}))
                error('Electrolyzer types defined for the modelling tool include only PEM and alkaline electrolysis.')
            end
            obj.type = string(lower(type));
        end
        
        
        function setElectrolyte(obj,electrolyte)
%           SETELECTROLYTE  Function for setting electrolyte for alkali electrolyzers.

            if strcmp(obj.type,"alkaline")
                if any(strcmpi(electrolyte,{'KOH','NaOH'}))
                    obj.electrolyte = electrolyte;
                else
                    error('Only KOH and NaOH defined as possible alkali electrolytes.')
                end
            else
                obj.electrolyte = "Polymer membrane";
            end
        end
        
        function potentialFunc = getPotential(obj,argin)
%            GETPOTENTIAL  Get the default func object for an user-defined potential.
%               Uses system data stored in the electrolyzerModel
%               object.
%               Input: argin -- Name of the potential as a string.
%               Output: potentialFunc -- default func object for the given
%                       potential.

            potentialName = strrep(string(lower(argin)), ' ', '');
            providedVariables = fieldnames(obj.Variables);
            switch string(lower(potentialName))
                case {"nernst","reversible","rev","opencircuit","ocv"}
                    if strcmpi(obj.type,'pem')
                        if all(ismember({'T','pCat','pAn'},providedVariables))
                            potentialFunc = nernst(obj.Variables.T,obj.Variables.pCat,obj.Variables.pAn,'type',obj.type);
                        else
                            error("To use Nernst equation with PEM the following variables, T (temperature in kelvin), pCat (cathode pressure in bar) and pAn (anode pressure in bar) have to be included in the electrolyzerModel Variables structure")
                        end
                    elseif strcmpi(obj.type,'alkaline')
                        if all(ismember({'T','p','m'},providedVariables))
                            potentialFunc = nernst(obj.Variables.T,obj.Variables.p,obj.Variables.m,'type',obj.type);
                        else
                            error("To use Nernst equation with alkaline the following variables, T (temperature in kelvin), p (csystem pressure in bar) and m (elektrolyte molality) have to be included in the electrolyzerModel Variables structure")
                        end
                    else
                        error("Electrolyzer type not recognised")
                    end
                case {"activation","act"}
                    potentialFunc = activation();
                case {"ohmic","ohm","resistive","res"}
                    potentialFunc = ohmic();
                case {"concentration","con","masstransfer"}
                    potentialFunc = concentration();
                otherwise
                    error("Potential component " + string(argin) + " not recognised.")
            end
        end
        
    end
    
end