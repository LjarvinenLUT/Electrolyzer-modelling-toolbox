classdef electrolyzerModel < model
    % ELECTROLYZERMODEL Contains all the necessary information about the electrolyzer model used and enables its simple use.
    %
    %   Herits properties and methods from MODEL class.
    %
    %   Instance of class ELECTROLYZERMODEL is used to create and store the
    %   function describing the UI curve of an electrolyzer, and store the
    %   experimental measurement results for its parameterisation.
    %   The class contains methods for simple fitting of the function to
    %   the measured data, as well as the calculation of the modelled
    %   voltage curve.
    %
    %   ELECTROLYZERMODEL Properties:
    %       electrolyte -- Chemical composition of the used electrolyte.
    %                       Significant to alkali electrolysis only.
    %       molarMassOfElectrolyte -- Molar mass of the electrolyte.
    %       shuntCurrentEnabled -- A boolean indicating the inclusion of
    %                               shunt current (current offset)
    %       type -- Type of the electrolyzer, PEM or alkaline
    %
    %   ELECTROLYZERMODEL Methods:
    %       addFunc -- Adds the expressions of the given potentials to
    %                       the modelFunc object.
    %       copy -- Creates a copy of the object whose properties are no
    %               longer linked to the parent.
    %       fitUI -- Performs a fit for the UI curve extracting values for
    %                   its parameters.
    %       mol2wtfrac -- Convert concentration from molality to weight
    %                       fractions.
    %       showUI -- Plots the UI curve for the model.
    %       shuntCurrent -- Enable or disable shunt current in the model.
    %       wtfrac2mol -- Convert concentration from weight fractions to
    %                       molality.
    %
    %   See also MODEL, FUNC 
    
    %% Properties
    properties (SetAccess = protected)
        type; % PEM or alkaline
        electrolyte; % Electrolyte
        molarMassOfElectrolyte; % Molar mass of the electrolyte, kg/mol
        shuntCurrentEnabled; % A boolean for inclusion of shunt current
    end
    
    
    %% Public methods
    methods
        
        function obj = electrolyzerModel(varargin)
            % ELECTROLYZERMODEL  Constructor function for the object.
            %  Optional inputs:
            %   type -- electrolyzer type.
            %              Default: PEM
            %   electrolyte -- Chemical composition of the electrolyte.
            %              Default: KOH (alkaline)
            %              Polymer membrane (PEM)
            %   nCells  -- Number of cells in series
            
            defaultType = "pem";
            defaultElectrolyte = "KOH";
            defaultNCells = 1;
            
            parser = inputParser;
            addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
            addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
            addParameter(parser,'nCells',defaultNCells,@(x) isnumeric(x)&&isscalar(x))
            
            parse(parser,varargin{:});
            
            parser.Results.type;
            
            % Create an empty model object
            obj@model;

            % Add number of cells in workspace
            obj.setInWorkspace(struct('Variables',...
                                      struct('nCells',...
                                             parser.Results.nCells)));
            
            setType(obj, parser.Results.type);
            setElectrolyte(obj, string(parser.Results.electrolyte));
        end
        

        function addFuncs(obj, varargin)
            % ADDFUNCS  Adds the given potential terms to modelFunc.
            %  Input of a list of string uses getPotential function to get
            %  the default func object. Alternatively the user can input a
            %  list of func object directly. Option 'rebuild' as the last
            %  parameter prevents addition of the components to
            %  funcStorage.
            %
            %   Overwrites MODEL.ADDFUNCS method
            %
            %  Examples:
            %
            %  obj.ADDFUNCS('nernst') adds the Nernst potential term to
            %   the modelFunc parameter using variables and
            %   electrolyzer type defined for the electrolyzerModel object.
            %
            %  obj.ADDFUNCS(func) adds the potential term defined by the
            %   input func object to modelFunc parameter.
            %
            %  obj.ADDFUNCS('nernst','activation',func1,func2) adds the
            %   all the given potential terms to the modelFunc
            %   parameter using the functionality meant for each type of
            %   input.
            %
            %  obj.ADDFUNCS(_,'names',nameCell) adds all the given
            %   potential terms to the modelFunc parameter and replaces
            %   their names in funcStorage with the ones given in a cell
            %   array after 'names' call.
            
            % Find entries that are of type char
            isCharEntries = cellfun(@ischar,varargin);
            indexOfCharEntries = find(isCharEntries);
            charEntries = varargin(isCharEntries);
            
            rebuildCall = ismember(charEntries,{'rebuild'});
            namesCall = ismember(charEntries,{'names'});
            
            % Use further the inputs that are not 'rebuild' nor 'names'
            otherInputs = true(size(varargin));
            
            % Check if rebuild mode was called
            if any(rebuildCall)
                rebuild = true;
                otherInputs(indexOfCharEntries(rebuildCall)) = false;
            else
                rebuild = false;
            end
            
            % Check if names were given for potentials
            if any(namesCall)
                if length(varargin) >= indexOfCharEntries(namesCall)+1 && (iscell(varargin{indexOfCharEntries(namesCall)+1})||isstring(varargin{indexOfCharEntries(namesCall)+1}))
                    names = varargin{indexOfCharEntries(namesCall)+1};
                    otherInputs(indexOfCharEntries(namesCall)+[0 1]) = false;
                else
                    error("Potential names have to be listed as a cell array after the 'names' parameter in the function call")
                end
            end
            
            inputs = varargin(otherInputs);
            
            % If only one input was given check if it is a cell or an array
            % and use its content if true.
            if length(inputs) == 1 && (iscell(inputs{1}) || (~iscell(inputs{1}) && length(inputs{1}) > 1))
                potentials = inputs{1};
            else
                potentials = inputs;
            end
            
            % Loop through all the potentials to be added
            for i = 1:length(potentials)
                if iscell(potentials)
                    potential = potentials{i};
                    if length(potentials(i))>1
                        if rebuild
                            obj.addFuncs(potential,'rebuild')
                        else
                            obj.addFuncs(potential)
                        end
                        continue;
                    end
                elseif isa(potentials,'func') % For func array inputs
                    potential = potentials(i);
                else
                    potential = potentials;
                end
                
                if isstring(potential) || ischar(potential) % String input
                    addedModelFunc = obj.getPotential(potential);
                    name = potential;
                elseif isa(potential,'func') % Func input
                    addedModelFunc = potential;
                    name = "unspecified";
                else
                    error("Potential to be added has to be a func object, or you have to specify with a string which potential term you want to add")
                end
                
                if any(namesCall) % If names were given in function call
                    try
                        name = names{i};
                    catch ME
                        if strcmp(ME.identifier,'MATLAB:badsubscript')
                            warning('Names were not specified for all the applied potentials. Please, check funcStorage property for the unspecified ones.')
                        else
                            rethrow(ME)
                        end
                    end
                end
                
                % Warn about multiply defined potential terms
                if isempty(addedModelFunc)
                    continue;
                elseif any(strcmp(addedModelFunc.equation,obj.funcStorage.equation)) && ~rebuild
                    warning("Same potential term included multiple times! For a list of added potential terms call the funcStorage property.")
                end
                
                
                obj.modelFunc = func.add(obj.modelFunc,addedModelFunc);                
                
                if ~rebuild % Supplement the funcStorage
                    StorageEntry = struct('name',name,'func',addedModelFunc,'equation',addedModelFunc.equation);
                    obj.funcStorage = [obj.funcStorage;struct2table(StorageEntry)];
                end
            end
        end


        function shuntCurrent(obj,state,printMsg)
            % SHUNTCURRENT Enables or disables shunt current.
            %   Input options:
            %       to enable shunt currents:
            %           true or "enable"
            %       to disable shunt currents:
            %           false or "disable"

            if ~exist("state","var") || isempty(state)
                % If state not provided, use the one already defined
                printMsg = false;
            elseif ismember(lower(string(state)),["enable" "enabled"]) || islogical(state) && state
                obj.shuntCurrentEnabled = true;
            elseif ismember(lower(string(state)),["disable" "disabled"]) || islogical(state) && ~state
                obj.shuntCurrentEnabled = false;
            else
                error("Unknown input "+string(state)+". Known inputs include true, false, 'enable' and 'disable'.")
            end

            if obj.shuntCurrentEnabled
                msg = "\nShunt currents enabled:\n"+...
                        "Model: constant shunt current\n";
            else
                msg = "\nShunt currents disabled\n";
            end

            % Print message if no value provided or the value is true
            if ~exist("printMsg","var") || islogical(printMsg) && printMsg
                fprintf(msg)
            end

            % Recreate the modelFunc to reset current
            obj.recreateModelFunc

            if obj.shuntCurrentEnabled
                SetWorkspace.Parameters.shuntCurrent = [];
                obj.setFitlims('shuntCurrent',{0,0,'min(x)'})
                obj.setInWorkspace(SetWorkspace)
                oldFuncHandle = func2str(obj.modelFunc.funcHandle);
                newFuncHandle = str2func(strrep(oldFuncHandle,...
                    'Workspace.Variables.current',...
                    '(Workspace.Variables.current-Workspace.Parameters.shuntCurrent)'));
                obj.modelFunc.setFuncHandle(newFuncHandle)
            else
                obj.removeFromWorkspace("shuntCurrent")
            end
        end

        

        
        function [fitParams,gof] = fitUI(obj,varargin)
            % FITUI Extracts the fit parameters for the electrolyzerModel.
            %  Calls external function fitUI.
            %
            %  obj.FITUI() performs UI curve fit using voltage and current
            %   data contained in the Variables structure of the
            %   electrolyzerModel object. Default method of particle
            %   swarm optimisation of the sum of square residuals is
            %   used for the fitting.
            %
            %  obj.FITUI(U,I) performs UI curve fit using voltage and
            %   current data provided as input parameters.
            %
            %  obj.FITUI(_,'method',m) performs UI curve fit using the
            %   method defined as input m.
            %   Options:
            %    'PS' -- Particle swarm optimisation
            %    'NLLSE' -- Non linear least squares error regression
            %  
            %  obj.FITUI(_,'weights',w) performs UI curve fit using the
            %   weighting method defined as input w.
            %   Options:
            %    'default' -- Weighs beginning and end of the curve
            %    'none' -- Doesn't add weights on the curve
            %               
            %  obj.FITUI(_,'plot',true) performs UI curve fit and plots the
            %   results.
            %
            %  Outputs:
            %   fitParams -- Table of fitting coefficient values and their
            %                   standard deviations
            %   gof -- Structure of goodnes-of-fit -parameters for the fit:
            %           - ssr -- Sum of Square Residuals
            %           - rmse -- Root Mean Square Error
            %           - rsquare -- R^2 value for the fit.
            %
            % See also FITUI
            
            defaultMethod = 'PS';
            defaultWeights = 'default';
            defaultPlot = false;
            defaultU = nan;
            defaultI = nan;
            
            parser = inputParser;
            addRequired(parser,'obj',@(x) isa(x,'electrolyzerModel'))
            addOptional(parser,'U',defaultU,@(x) isnumeric(x))
            addOptional(parser,'I',defaultI,@(x) isnumeric(x))
            addParameter(parser,'method',defaultMethod,@(x) ischar(x)||isstring(x)) % Fitting method to be used
            addParameter(parser,'weights',defaultWeights,@(x) ischar(x)||isstring(x))
            addParameter(parser,'plot',defaultPlot,@(x) islogical(x))
            
            parse(parser,obj,varargin{:});
            
            U = parser.Results.U;
            I = parser.Results.I;
            method = upper(string(parser.Results.method));
            weightsMethod = lower(string(parser.Results.weights));
            usePlotting = parser.Results.plot;

            % Check that if multiple models are provided they are in a 1xn
            % vector
            if min(size(obj)) > 1
                error("Model objects should be provided as an 1xn vector")
            end

            nModels = numel(obj); % Number of models to be fit
            fitParams = cell(size(obj));
            gof = cell(size(obj));

            for i = 1:nModels
                % Apply shunt current if enabled
                obj(i).shuntCurrent;

                if isnan(I) % Current not provided as an input
                    if any(ismember('current',fieldnames(obj(i).modelFunc.Workspace.Variables)))
                        I = obj(i).modelFunc.Variables.current;
                    else
                        error("Current not defined. Defined it either in the Variables structure of the electrolyzerModel object. Alternatively, if voltage has been provided as the second parameter for fitUI, current can be provided as the third parameter.")
                    end
                    if isnan(U) % Voltage not provided as an input
                        if any(ismember('voltage',fieldnames(obj(i).modelFunc.Workspace.Variables)))
                            U = obj(i).modelFunc.Variables.voltage;
                        else
                            error("Voltage not defined. Defined it either in the Variables structure of the electrolyzerModel object. Alternatively, if voltage has been provided as the second parameter for fitUI, current can be provided as the third parameter.")
                        end
                    end
                end

                [fitParams{i},gof{i}] = obj.performFit(obj(i).modelFunc,U,I,'method',method,'weights',weightsMethod);

                obj(i).PlottingCurves = struct('currentMeasured',I,'voltageMeasured',U,'method',method);

                obj(i).synchronizeFuncStorage;

                if usePlotting
                    obj(i).showUI;
                end
            end

            if nModels == 1
                fitParams = fitParams{1};
                gof = gof{1};
            end
            
        end
        
        function result = calculate(obj,varargin)
            % CALCULATE Calculates voltage from the UI curve.
            %   Calls the calculate method of the func object
            %   modelFunc.
            %   Variables input as name-value pairs are used for the
            %   calculation and any variable required by the function that 
            %   is not input to CALCULATE is looked for from the Workspace 
            %   of the modelFunc object.
            %
            %   See also FUNC.CALCULATE
            
            if min(size(obj)) > 1
                error("Model objects should be provided as a vector array")
            end

            nModels = numel(obj); % Number of models to be calculated
            result = cell(size(obj));
            for i = 1:nModels
                tempResult = obj(i).modelFunc.calculate(varargin{:});
                if any(~isreal(tempResult))
                    warning("Complex voltage values detected. This may indicate that the values of the given current vector exceed the fit limitting current density 'j_lim'. Imaginary parts ignored.")
                end
                result{i} = real(tempResult);
            end
            
            if nModels == 1
                result = result{1};
            end

        end
        
        
        function showUI(obj,varargin)
            % SHOWUI Creates a figure and plots the UI curve on it.
            %   Uses data stored in PlottingCurves property by the fitUI
            %   method.
            %
            %   Optional input: 
            %       'label' - Text to be showed as the name of the
            %                   created uifigure window. (Name-value pair)
            %       'separate' - Show separate potential components as a
            %                   stacked plot
            %
            %   Examples:
            %
            %   electrolyzerModel.SHOWUI();
            %
            %   electrolyzerModel.SHOWUI(_,'label','message') -- Create a
            %       figure with 'message' as the figure name.
            %
            %   electrolyzerModel.SHOWUI(_,'separate') -- Create a UI plot
            %       with the voltage components separated 
            %       TODO: use patch (https://se.mathworks.com/help/matlab/ref/patch.html)
            %

            textInput = isstring(varargin)||ischar(varargin);
            if nargin == 2 && (isstring(varargin{1})||ischar(varargin{1}))
                figureMsg = varargin{1};
            else
                figureMsg = 'Automatically created plot of the UI curve and its parameters';
            end
            
            if isempty(obj.PlottingCurves)
                warning('No UI curve fit performed, therefore no UI curve to show.')
            else
                current = obj.PlottingCurves.currentMeasured;
                voltage = obj.PlottingCurves.voltageMeasured;
                
                % Calculate dense data vectors for the fitted UI curve
                fullFitCurrent = min(current):0.001:max(current);
                fullFitVoltage = mean(obj.calculate('current',fullFitCurrent),1);
                
                % Take samples from the dense data vectors
                N = 100; % Number of evenly taken voltage samples
                voltageSamples = linspace(min(fullFitVoltage),max(fullFitVoltage),N)';
                for ii = 1:N
                    Udif = abs(fullFitVoltage-voltageSamples(ii));
                    [~,ind] = min(Udif);
                    fitCurrent(ii,1) = fullFitCurrent(ind); % Final sampled current vector
                    fitVoltage(ii,1) = fullFitVoltage(ind); % Final sampled voltage vector
                end
                
                obj.PlottingCurves.currentFit = fitCurrent;
                obj.PlottingCurves.voltageFit = fitVoltage;
                
                % Build a cell array for errorbar input
                if length(voltage(1,:)) == 1 % No voltage error given
                    voltageplot = {voltage,[],[]};
                elseif length(voltage(1,:)) == 2 
                    voltageplot = {voltage(:,1),voltage(:,2),voltage(:,2)};
                elseif length(voltage(1,:)) > 2
                    msg = ['Too wide voltage data matrix! Provide the data'...
                        ' either as a column vector or as a nx2 matrix'...
                        ' with standard deviation in the second column'];
                    error(msg)
                else 
                    error('No voltage to plot!')
                end
                
                if length(current(1,:)) == 1 % No current error given
                    cvplot = [current,voltageplot,{[],[]}];
                elseif length(current(1,:)) >= 2
                    cvplot = [current(:,1),voltageplot,current(:,2),current(:,2)];
                elseif length(current(1,:)) > 2
                    msg = ['Too wide current data matrix! Provide the data'...
                        ' either as a column vector or as a nx2 matrix'...
                        ' with standard deviation in the second column'];
                    error(msg)
                else
                    error('No current to plot!')
                end
%               
                % Figure parameters
                dummyFig = figure;
                set(dummyFig,'Units','pixels')
                initFigPosition = get(dummyFig,'Position'); % Initial figure dimensions [left bottom width height]
                close(dummyFig)
                f = uifigure('Name',figureMsg,'HandleVisibility', 'on','Position',initFigPosition);
                
                % Resizing the figure
                newHeight = initFigPosition(4)*1.5;
                newBottom = max(0,f.Position(2)+f.Position(4)-newHeight);
                f.Position([2,4]) = [newBottom newHeight]; % Determine new dimensions
                
                figMeas = f.Position;
                
                format shortG % Format the table to show small numbers as the powers of ten
                
                % Add the UI curve plot
                ax = uiaxes(f);
                ax.Position = [figMeas(3)*0.025 figMeas(4)*0.325 figMeas(3)*0.95 figMeas(4)*0.65];
                errorbar(ax,cvplot{:},'o')
                hold(ax,'on');
                plot(ax,fitCurrent,fitVoltage)
                hold(ax,'off');
                textInterpreter = get(groot,'defaultTextInterpreter');
                if strcmp(textInterpreter,'latex')
                    ax.XLabel.String = "Current density (A/cm$^2$)";
                else
                    ax.XLabel.String = "Current density (A/cm^2)";
                end
                ax.YLabel.String = "Voltage (V)";
                legend(ax,"Data", "Fit", "Location", "Best")
                
                % Add label about the used method
                txa = uitextarea(f);
                txa.Position = [figMeas(3)*0.025 figMeas(4)*0.275 figMeas(3)*0.95 figMeas(4)*0.05];
                if strcmpi(obj.PlottingCurves.method,"PS")
                    methodStr = "Particle swarm minimization on the sum of squared residuals";
                elseif strcmpi(obj.PlottingCurves.method,"NLLSE")
                    methodStr = "Non-linear least squares error";
                end 
                txa.Value = "Fit method: " + methodStr;
                
                % Add fit parameters to the plot
                paramTable = obj.getParams;
                data = table2cell(paramTable);
                uit = uitable(f);
                uit.Position = [figMeas(3)*0.025 figMeas(4)*0.025 figMeas(3)*0.95 figMeas(4)*0.225];
                uit.Data = data;
                uit.ColumnName = {'Value','STD'};
                uit.RowName = paramTable.Properties.RowNames;
                uit.ColumnFormat = {'shortG','shortG'}; % Format to scientific notation if better
            end
        end


        
        function report(obj)
            % REPORT Displays a report of all the properties of the object
            %
            % See also MODEL.REPORT, MODEL.VIEWWORKSPACE
            msg1 = ['\nElectrolyzer model properties:\n'...
                    '\n type: %s\n' ...
                    '\n electrolyte: %s\n' ...
                    '\n molarMassOfElectrolyte: %6.4f kg/mol\n'];
            fprintf(msg1,obj.type,obj.electrolyte,obj.molarMassOfElectrolyte)
            report@model(obj) %
        end
        
        
        function childObj = copy(obj)
            % COPY  Creates a full copy of the object with its own handle.
            %   Properties of the child object are no longer related to
            %   those of the parent.
            childObj = electrolyzerModel('type',obj.type,'electrolyte',obj.electrolyte);
            if ~func.isEmpty(obj.modelFunc)
                childObj.addFuncs(obj.funcStorage.func.copy,'names',obj.funcStorage.name);
            end
            childObj.setInWorkspace(obj.modelFunc.Workspace);
            childObj.shuntCurrent(obj.shuntCurrentEnabled,false);
        end
        
        function molality = wtfrac2molal(obj,wtfrac)
            % WTFRAC2MOLAL Convert concentration value in
            %   weight fraction (mass of solute/mass of solution) to
            %   molality (moles of solute/mass of solvent in kg) using the
            %   molar mass value contained in the instance.
            %
            %   See also WTFRAC2MOL
            molality = wtfrac2molal(wtfrac,obj.molarMassOfElectrolyte);
        end
        
        function wtfrac = molal2wtfrac(obj,molality)
            % MOLAL2WTFRAC Convert concentration value in
            %   molality (moles of solute/mass of solvent in kg) to
            %   weight fraction (mass of solute/mass of solution) using the
            %   molar mass value contained in the instance.
            %
            %   See also MOL2WTFRAC
            wtfrac = molal2wtfrac(molality,obj.molarMassOfElectrolyte);
        end
        
    end
    
    %% Private methods
    methods (Access = private)
        
        function setType(obj,type)
            % SETTYPE  Sets the electrolyzer type.
            if strcmpi(type,'alkali')
                type = 'alkaline';
            elseif ~any(strcmpi(type,{'pem','alkaline'}))
                error('Electrolyzer types defined for the modelling tool include only PEM and alkaline electrolysis.')
            end
            obj.type = string(lower(type));
        end
        
        
        function setElectrolyte(obj,electrolyte)
            % SETELECTROLYTE Sets the electrolyte for alkali electrolyzers.
            if strcmp(obj.type,"alkaline")
                switch electrolyte
                    case 'KOH'
                        obj.molarMassOfElectrolyte = 39.0983 + 15.9994 + 1.0079;
                        obj.setInWorkspace(struct(...
                                'Variables',struct(...
                                    'electrolyte',1,...
                                    'molarMassOfElectrolyte',obj.molarMassOfElectrolyte,...
                                    'molality',[],...
                                    'Molarity',[],...
                                    'wtfrac',[],...
                                    'T',[])));
                    case 'NaOH'
                        obj.molarMassOfElectrolyte = 22.9898 + 15.9994 + 1.0079;
                        obj.setInWorkspace(struct(...
                                'Variables',struct(...
                                    'electrolyte',2,...
                                    'molarMassOfElectrolyte',obj.molarMassOfElectrolyte,...
                                    'molality',[],...
                                    'Molarity',[],...
                                    'wtfrac',[],...
                                    'T',[])));
                    otherwise
                        error('Only KOH and NaOH defined as possible alkali electrolytes.')
                end
                obj.setInWorkspace(struct(...
                                'Dependencies',struct(...
                                    'Molarity',"if changed.Molarity;"+...
                                                "Workspace.Variables.molality = molar2molal(Workspace.Variables.Molarity,Workspace.Variables.T,'"+string(electrolyte)+"');"+...
                                                "Workspace.Variables.wtfrac = molal2wtfrac(Workspace.Variables.molality,Workspace.Variables.molarMassOfElectrolyte);"+...
                                                "end;",...
                                    'molality',"if changed.molality;"+...
                                                "Workspace.Variables.Molarity = molal2molar(Workspace.Variables.molality,Workspace.Variables.T,'"+string(electrolyte)+"');"+...
                                                "Workspace.Variables.wtfrac = molal2wtfrac(Workspace.Variables.molality,Workspace.Variables.molarMassOfElectrolyte);"+...
                                                "end;",...
                                    'wtfrac',"if changed.wtfrac;"+...
                                                "Workspace.Variables.molality = wtfrac2molal(Workspace.Variables.wtfrac,Workspace.Variables.molarMassOfElectrolyte);"+...
                                                "Workspace.Variables.Molarity = molal2molar(Workspace.Variables.molality,Workspace.Variables.T,'"+string(electrolyte)+"');"+...
                                                "end;",...
                                    'concentration_temperature',"if changed.T;"+...
                                                "Workspace.Variables.Molarity = molal2molar(Workspace.Variables.molality,Workspace.Variables.T,'"+string(electrolyte)+"');"+...
                                                "Workspace.Variables.wtfrac = molal2wtfrac(Workspace.Variables.molality,Workspace.Variables.molarMassOfElectrolyte);"+...
                                                "end;")));
                obj.electrolyte = electrolyte;
            else
                obj.electrolyte = "Polymer membrane";
            end
        end
        
        function potentialFunc = getPotential(obj,argin)
            % GETPOTENTIAL  Get the default func object for an user-defined potential.
            %  Uses system data stored in the electrolyzerModel object.
            %  Input: argin -- Name of the potential as a string.
            %  Output: potentialFunc -- default func object for the given
            %                           potential.
            potentialName = strrep(string(lower(argin)), ' ', '');
%             providedVariables = fieldnames(obj.potentialFunc.Workspace.Variables);
            switch string(lower(potentialName))
                case {"nernst","reversible","rev","opencircuit","ocv"}
%                     if strcmpi(obj.type,'pem')
%                         if all(ismember({'T','pCat','pAn'},providedVariables))
%                             potentialFunc = nernst(obj.type);
%                         else
%                             error("To use Nernst equation with PEM the following variables, T (temperature in kelvin), pCat (cathode pressure in bara) and pAn (anode pressure in bara) have to be included in the electrolyzerModel Variables structure")
%                         end
%                     elseif strcmpi(obj.type,'alkaline')
%                         if all(ismember({'T','ps','m'},providedVariables))
%                             potentialFunc = nernst(obj.type);
%                         else
%                             error("To use Nernst equation with alkaline the following variables, T (temperature in kelvin), ps (system pressure in bara) and m (electrolyte molality) have to be included in the electrolyzerModel Variables structure. The variables have to have the exact same naming as shown in the previous sentence.")
%                         end
%                     else
%                         error("Electrolyzer type not recognised")
%                     end
                    potentialFunc = nernst(obj.type);
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