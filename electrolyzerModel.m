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
    %       funcStorage -- A table that stores the separate func objects
    %                       with their most significant information
    %                       visible. Storage enables the removal of
    %                       unnecessary potential components and automatic
    %                       check for doubly included ones.
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
    %       addPotentials -- Adds the expressions of the given potentials to
    %                       the potentialFunc object.
    %       calculate -- Calculates voltage from the UI curve for given 
	%                       currents.
    %       clearOverpotentials -- Clears the potentialFunc property.
    %       copy -- Creates a copy of the object whose properties are no
    %               longer linked to the parent.
    %       fitUI -- Performs a fit for the UI curve extracting values for
    %                   its coefficients.
    %       getCoefficients -- Outputs the coefficient structure from
    %                           potentialFunc Workspace.
    %       removeVariables -- Remove variables from the Variables
    %                           structure.
    %       replaceParams --  A method for replacing parameters in the 
    %                           potentialFunc Workspace structure.
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
       funcStorage; % A table that stores the separate func objects
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
            obj.Variables = struct();
            
            % Create funcStorage
            variableNamesTypes = [["name", "string"];...
                                    ["func", "func"];...
                                    ["equation", "string"]];
            obj.funcStorage = table('Size',[0,size(variableNamesTypes,1)],... 
                                    'VariableNames', variableNamesTypes(:,1),...
                                    'VariableTypes', variableNamesTypes(:,2));
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
        
        
        function replaceParams(obj,varargin)
%           REPLACEPARAMS  A method for replacing parameters in the potentialFunc Workspace structure.
%               Parameters should be provided as a name-value pair or 
%               directly as a structure.
            if isempty(varargin{1})
                return;
            elseif length(varargin) == 1 && isstruct(varargin{1})
                paramsToReplace = varargin{1};
            elseif mod(nargin,2)
                for i = 1:2:length(varargin)
                    paramsToReplace.(varargin{i}) = varargin{i+1};
                end
            else
                error('Parameters to be replaced have to be either set as a single structure or as name-value pairs')
            end
            
            obj.potentialFunc.Workspace = addValuesToStruct(obj.potentialFunc.Workspace,paramsToReplace);
        end
        
        
        
        function addPotentials(obj, varargin)
%           ADDPOTENTIALS  Adds the given potential terms to the total potential func object. 
%               Input of a list of string uses getPotential function to get the 
%               default func object. Alternatively the user can input a 
%               list of func object directly. Option 'rebuild' as the last
%               parameter 
%
%           Examples:
%
%           obj.ADDPOTENTIALS('nernst') adds the Nernst potential term to
%               the potentialFunc parameter using variables and 
%               electrolyzer type defined for the electrolyzerModel object.
%
%           obj.ADDPOTENTIALS(func) adds the potential term defined by the
%               input func object to potentialFunc parameter.
%
%           obj.ADDPOTENTIALS('nernst','activation',func1,func2) adds the
%               all the given potential terms to the potentialFunc
%               parameter using the functionality meant for each type of
%               input.

            if (isstring(varargin{end})||ischar(varargin{end}))&&strcmp(varargin{end},'rebuild')
                rebuild = true;
                inputs = varargin{1:end-1};
            else
                rebuild = false;
                inputs = varargin;
            end
            if length(inputs) == 1 && iscell(inputs{1})
                potentials = inputs{1};
            else
                potentials = inputs;
            end
            
            for i = 1:length(potentials)
                if iscell(potentials)
                    potential = potentials{i};
                elseif isa('potentials','func') % For func array inputs
                    potential = potentials(i);
                end
                
                if isstring(potential) || ischar(potential)
                    addedPotentialFunc = obj.getPotential(potential);
                    name = potential;
                elseif isa(potential,'func')
                    addedPotentialFunc = potential;
                    name = "unspecified";
                else
                    error("Potential to be added has to be a func object, or you have to specify with a string which potential term you want to add")
                end
                
                if ~isa(obj.potentialFunc,'func') % If no previous potential function is assigned
                    obj.potentialFunc = addedPotentialFunc;
                else
                    obj.potentialFunc = func.add(obj.potentialFunc,addedPotentialFunc);
                end
                
                if ~rebuild
                    StorageEntry = struct('name',name,'func',addedPotentialFunc,'equation',addedPotentialFunc.equation);
                    obj.funcStorage = [obj.funcStorage;struct2table(StorageEntry)];
                end
            end
        end
        
        function removePotentials(obj,varargin)
%           REMOVEPOTENTIALS  Removes the given potential terms from the total potential func object. 
%               Input can be either in the form of string for removing
%               potentials of certain name, or numeric for removing
%               potentials in certain indeces. Multiple potentials can be
%               removed by listing them separately or as a cell array.
%
%           Examples:
%
%           obj.REMOVEPOTENTIALS('nernst') removes the Nernst potential term
%               from the potentialFunc parameter.
%
%           obj.REMOVEPOTENTIALS(index) removes the potential term with the
%               given index in funcStorage.
%
%           obj.REMOVEPOTENTIALS('nernst','activation',ind1,ind2) removes
%               all the given potential terms from the potentialFunc
%               parameter using the functionality meant for each type of
%               input.

            if length(varargin) == 1 && iscell(varargin{1})
                potentials = varargin{1};
            else
                potentials = varargin;
            end
            
            % Remove the given potentials from the storage
            for i = 1:length(potentials)
                if isstring(potentials{i}) || ischar(potentials{i})
                    ind = find(strcmp(obj.funcStorage.name,potentials{i}));
                    obj.funcStorage(ind,:) = [];
                elseif isnumeric(potentials{i})
                    obj.funcStorage(potentials{i},:) = [];
                else
                    error("Potential to be removed not recogniced. Input either as a string for the name or number for its index in funcStorage.")
                end
            end
            
            % Recreate potentialFunc from the storage
            obj.clearPotentials;
            funcs = obj.funcStorage.func;
            obj.addPotentials(funcs,'rebuild');
            
        end
        
         
        function clearPotentials(obj)
            % CLEARPOTENTIALS Clears the potential function
            obj.potentialFunc = "Potential function not defined";
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
        
        function Coefficients = getCoefficients(obj)
            % GETCOEFFICIENTS Outputs a structure containing fit coefficients
            Coefficients = obj.potentialFunc.Workspace.Coefficients;
        end
        
        function varargout = viewWorkspace(obj)
            % VIEWWORKSPACE Outputs the workspace of potentialFunc in a human-readable table
            report = obj.potentialFunc.viewWorkspace;
            disp(report)
            if nargout == 1
                varargout{1} = report;
            elseif nargout > 1
                error("Too many output arguments.")
            end
        end
        
        function childObj = copy(obj)
            % COPY  Creates a full copy of the object with its own handle.
            %   Properties of the child object are no longer related to
            %   those of the parent.
            childObj = electrolyzerModel('type',obj.type,'electrolyte',obj.electrolyte);
            if isa(obj.potentialFunc,'func')
                childObj.addPotentials(obj.potentialFunc.copy);
            end
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