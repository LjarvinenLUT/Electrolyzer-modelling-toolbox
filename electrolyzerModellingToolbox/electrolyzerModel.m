classdef electrolyzerModel < handle
    % ELECTROLYZERMODEL Contains all the necessary information about the electrolyzer model used and enables its simple use.
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
    %       funcStorage -- A table that stores the separate func objects
    %                       with their most significant information
    %                       visible. Storage enables the removal of
    %                       unnecessary potential components and automatic
    %                       check for doubly included ones.
    %       molarMassOfElectrolyte -- Molar mass of the electrolyte.
    %       PlottingCurves -- A structure containing the measurements used
    %                           for curve fitting as well as the calculated
    %                           curve from the fit results
    %       potentialFunc -- A func object containing the UI curve
    %                           determining equation, values of its
    %                           variables, constants and parameters, and
    %                           functionality for UI fitting and
    %                           calculation.
    %       shuntCurrentEnabled -- A boolean indicating the inclusion of
    %                               shunt current (current offset)
    %       type -- Type of the electrolyzer, PEM or alkaline
    %
    %   ELECTROLYZERMODEL Methods:
    %       addPotentials -- Adds the expressions of the given potentials to
    %                       the potentialFunc object.
    %       calculate -- Calculates voltage from the UI curve for given
    %                       currents.
    %       clearPotentials -- Clears the potentialFunc property.
    %       copy -- Creates a copy of the object whose properties are no
    %               longer linked to the parent.
    %       fitUI -- Performs a fit for the UI curve extracting values for
    %                   its parameters.
    %       getParams -- Outputs the parameters from
    %                           potentialFunc Workspace in a table.
    %       mol2wtfrac -- Convert concentration from molality to weight
    %                       fractions.
    %       removePotentials -- Remove given potential terms from the
    %                           potentialFunc object and the func Storage.
    %       removeFromWorkspace -- Remove values from the potentialFunc
    %                           Workspace structure.
    %       replaceInWorkspace --  A method for replacing values in
    %                               the potentialFunc Workspace structure.
    %       report -- Print a report of all the properties of the modelling
    %                   object.
    %       showUI -- Plots the UI curve for the model.
    %       setFitlims -- Set fit limits for the model parameters.
    %       setInWorkspace -- Set values in the Workspace structure of
    %                           the potentialFunc object.
    %       shuntCurrent -- Enable or disable shunt current in the model.
    %       viewWorkspace -- Outputs the workspace of potentialFunc in a
    %                           human-readable table.
    %       wtfrac2molal -- Convert concentration from weight fractions to
    %                       molality.
    %
    %   See also FUNC
    
    %% Properties
    properties (SetAccess = protected)
        type; % PEM or alkaline
        electrolyte; % Electrolyte
        molarMassOfElectrolyte; % Molar mass of the electrolyte, kg/mol
        potentialFunc; % A func object for combined overpotential function
        funcStorage; % A table that stores the separate func objects
        PlottingCurves; % A structure with the measurements and fit curve
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
            
            defaultType = "pem";
            defaultElectrolyte = "KOH";
            
            parser = inputParser;
            addParameter(parser,'type',defaultType,@(x) ischar(x)||isstring(x));
            addParameter(parser,'electrolyte',defaultElectrolyte,@(x) ischar(x)||isstring(x))
            
            parse(parser,varargin{:});
            
            parser.Results.type;
            
            % Create an empty func object
            obj.potentialFunc = func.createEmpty;
            
            setType(obj, parser.Results.type);
            setElectrolyte(obj, string(parser.Results.electrolyte));
            
            % Create funcStorage
            variableNamesTypes = [["name", "string"];...
                ["func", "func"];...
                ["equation", "string"]];
            obj.funcStorage = table('Size',[0,size(variableNamesTypes,1)],...
                'VariableNames', variableNamesTypes(:,1),...
                'VariableTypes', variableNamesTypes(:,2));
            % Create empty structure for plotting curves
            obj.PlottingCurves = struct([]);
        end
        
        function setInWorkspace(obj,SetWorkspace)
            % SETINWORKSPACE  Sets values in the Workspace of potentialFunc.
            %  Parameters to be set should be provided as a Workspace
            %  structure. To check compatibility, use func.isWorkspace.
            %
            %  See also FUNC.SETINWORKSPACE, FUNC.ISWORKSPACE
            obj.potentialFunc.setInWorkspace(SetWorkspace)
            obj.synchronizeFuncStorage;
        end
        
        function removeFromWorkspace(obj,varargin)
            % REMOVEFROMWORKSPACE Removes parameters from the Workspace structure.
            %  Parameters to be removed have to be provided either as a list or
            %  a cell of strings/character vectors.
            %
            %  See also FUNC.REMOVEFROMWORKSPACE
            if length(varargin) == 1 && iscell(varargin{1})
                valuesToRemove = varargin{1};
            else
                valuesToRemove = varargin;
            end
            
            obj.potentialFunc.removeFromWorkspace(valuesToRemove)
            obj.synchronizeFuncStorage;
        end
        
        
        function replaceInWorkspace(obj,varargin)
            % REPLACEINWORKSPACE Replaces parameters in the Workspace structure.
            %  Parameters should be provided as a name-value pair or directly
            %  as a structure.
            %
            % See also FUNC.REPLACEINWORKSPACE
            obj.potentialFunc.replaceInWorkspace(varargin{:});
            obj.synchronizeFuncStorage;
        end
        
        
        
        function addPotentials(obj, varargin)
            % ADDPOTENTIALS  Adds the given potential terms to potentialFunc.
            %  Input of a list of string uses getPotential function to get
            %  the default func object. Alternatively the user can input a
            %  list of func object directly. Option 'rebuild' as the last
            %  parameter prevents addition of the components to
            %  funcStorage.
            %
            %  Examples:
            %
            %  obj.ADDPOTENTIALS('nernst') adds the Nernst potential term to
            %   the potentialFunc parameter using variables and
            %   electrolyzer type defined for the electrolyzerModel object.
            %
            %  obj.ADDPOTENTIALS(func) adds the potential term defined by the
            %   input func object to potentialFunc parameter.
            %
            %  obj.ADDPOTENTIALS('nernst','activation',func1,func2) adds the
            %   all the given potential terms to the potentialFunc
            %   parameter using the functionality meant for each type of
            %   input.
            %
            %  obj.ADDPOTENTIALS(_,'names',nameCell) adds all the given
            %   potential terms to the potentialFunc parameter and replaces
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
                            obj.addPotentials(potential,'rebuild')
                        else
                            obj.addPotentials(potential)
                        end
                        continue;
                    end
                elseif isa(potentials,'func') % For func array inputs
                    potential = potentials(i);
                else
                    potential = potentials;
                end
                
                if isstring(potential) || ischar(potential) % String input
                    addedPotentialFunc = obj.getPotential(potential);
                    name = potential;
                elseif isa(potential,'func') % Func input
                    addedPotentialFunc = potential;
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
                if any(strcmp(addedPotentialFunc.equation,obj.funcStorage.equation)) && ~rebuild
                    warning("Same potential term included multiple times! For a list of added potential terms call the funcStorage property.")
                end
                
                
                obj.potentialFunc = func.add(obj.potentialFunc,addedPotentialFunc);                
                
                if ~rebuild % Supplement the funcStorage
                    StorageEntry = struct('name',name,'func',addedPotentialFunc,'equation',addedPotentialFunc.equation);
                    obj.funcStorage = [obj.funcStorage;struct2table(StorageEntry)];
                end
            end
        end
        
        
        function removePotentials(obj,varargin)
            % REMOVEPOTENTIALS  Removes potential terms from potentialFunc.
            %  Input can be either in the form of string for removing
            %  potentials of certain name, or numeric for removing
            %  potentials in certain indices in the funcStorage. Multiple
            %  potentials can be removed by listing them separately or as
            %  a cell array.
            %
            %  Refer to property funcStorage to find the available names
            %  and indeces.
            %
            %  Examples:
            %
            %  obj.REMOVEPOTENTIALS('nernst') removes the Nernst potential term
            %   from the potentialFunc parameter.
            %
            %  obj.REMOVEPOTENTIALS(index) removes the potential term with the
            %   given index in funcStorage.
            %
            %  obj.REMOVEPOTENTIALS('nernst','activation',ind1,ind2) removes
            %   all the given potential terms from the potentialFunc
            %   parameter using the functionality meant for each type of
            %   input.
            
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
            obj.recreatePotentialFunc

        end
        
        
        function clearPotentials(obj,varargin)
            % CLEARPOTENTIALS Clears the potential function.
            %   Replaces it with an empty FUNC object.
            %
            % See also FUNC.CREATEEMPTY
            
            obj.potentialFunc = obj.potentialFunc.copy('empty');
            if nargin == 2 && strcmp(varargin{1},'rebuild')
                return;
            else
                obj.funcStorage(1:height(obj.funcStorage),:) = [];
            end
        end
        
        
        
        function setFitlims(obj,varargin)
            % SETFITLIMS Sets the potentialFunc property fitLims.
            %   Set fitting limit values for the model parameters.
            %   Fitlims structure contains one field for each coefficient,
            %   named identical to the coefficient. Each field contains a
            %   1x3 cell array with values {low,start,high} for lower
            %   limit, starting point and higher limit, respectively.
            %
            %   Recognized input is a Structure with the right format or
            %   name value pairs with the coefficient name followed by the
            %   limit matrix.
            %
            %   Limits can be either as:
            %       - numeric scalars or
            %       - MATLAB equations in string form each outputing a
            %       numeric scalar based on the dependent variable values.
            %       Dependent variable has to be signed with 'x' in the
            %       equations.
            %
            %   See also FUNC.SETFITLIMS
            
            obj.potentialFunc.setFitlims(varargin{:});
        end



        function shuntCurrent(obj,state)
            % SHUNTCURRENT Enables or disables shunt current.
            %   Input options:
            %       to enable shunt currents:
            %           true or "enable"
            %       to disable shunt currents:
            %           false or "disable"

            if ~exist("state","var")
                % If state not provided, use the one already defined
            elseif ismember(lower(string(state)),["enable" "enabled"]) || islogical(state) && state
                obj.shuntCurrentEnabled = true;
                fprintf('\nShunt currents enabled:\n')
                fprintf('Model: constant shunt current\n')
            elseif ismember(lower(string(state)),["disable" "disabled"]) || islogical(state) && ~state
                obj.shuntCurrentEnabled = false;
                fprintf('\nShunt currents disabled\n')
            else
                error("Unknown input "+string(state)+". Known inputs include true, false, 'enable' and 'disable'.")
            end

            % Recreate the potentialFunc to reset current
            obj.recreatePotentialFunc

            if obj.shuntCurrentEnabled
                SetWorkspace.Parameters.shuntCurrent = [];
                obj.setFitlims('shuntCurrent',{0,0,'min(x)'})
                obj.setInWorkspace(SetWorkspace)
                oldFuncHandle = func2str(obj.potentialFunc.funcHandle);
                newFuncHandle = str2func(strrep(oldFuncHandle,...
                    'Workspace.Variables.current',...
                    '(Workspace.Variables.current-Workspace.Parameters.shuntCurrent)'));
                obj.potentialFunc.setFuncHandle(newFuncHandle)
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
                    if any(ismember('current',fieldnames(obj(i).potentialFunc.Workspace.Variables)))
                        I = obj(i).potentialFunc.Variables.current;
                    else
                        error("Current not defined. Defined it either in the Variables structure of the electrolyzerModel object. Alternatively, if voltage has been provided as the second parameter for fitUI, current can be provided as the third parameter.")
                    end
                    if isnan(U) % Voltage not provided as an input
                        if any(ismember('voltage',fieldnames(obj(i).potentialFunc.Workspace.Variables)))
                            U = obj(i).potentialFunc.Variables.voltage;
                        else
                            error("Voltage not defined. Defined it either in the Variables structure of the electrolyzerModel object. Alternatively, if voltage has been provided as the second parameter for fitUI, current can be provided as the third parameter.")
                        end
                    end
                end

                [fitParams{i},gof{i}] = fitUI(obj(i).potentialFunc,U,I,'method',method,'weights',weightsMethod);

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
            %   potentialFunc.
            %   Variables input as name-value pairs are used for the
            %   calculation and any variable required by the function that 
            %   is not input to CALCULATE is looked for from the Workspace 
            %   of the potentialFunc object.
            %
            %   See also FUNC.CALCULATE
            
            if min(size(obj)) > 1
                error("Model objects should be provided as a vector array")
            end

            nModels = numel(obj); % Number of models to be calculated
            result = cell(size(obj));
            for i = 1:nModels
                tempResult = obj(i).potentialFunc.calculate(varargin{:});
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
            %   Optional input: Text to be showed as the name of the
            %       created uifigure window.
            
            % 
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
        
        % TODO: should be made more robust. Checking missing variables for
        % the calculation should be automatic.
        function out =  showOverpotentials(obj, varargin)
            % SHOWOVERPOTENTIALS plots individual overpotentials
            %  obj.SHOWOVERPOTENTIALS() plots the overpotentials for
            %  current range of the measurements
            %  obj.SHOWOVERPOTENTIALS(current) plots the overpotentials for
            %  current range specified by the array
            if nargin() == 1
                current = obj.PlottingCurves.currentMeasured;
            else
                current = cell2mat(varargin(1));
            end
            out = {};
            funcCount = length(obj.funcStorage.func);
            figure()
            hold on;
            for i = 1:funcCount
                % Dont plot OCV
                if obj.funcStorage.name(i) ~= "nernst" && obj.funcStorage.name(i) ~= "ocv"
                    plot(current, obj.funcStorage.func(i).calculate('current', current));
                    out(end + 1) = {obj.funcStorage.func(i).calculate('current', current)};
                end
            end
            hold off;
            % Exclude names ocv and nernst from the legend text
            filteredNames = obj.funcStorage.name(find(obj.funcStorage.name ~= "nernst" & ...
                obj.funcStorage.name ~= "ocv"));
            legend(filteredNames)
            xlabel("Current density (A/cm$^2$)")
            ylabel("Overpotential (V)")
        end
        
        function coeffTable = getParams(obj)
            % GETPARAMS Outputs a table containing fit parameters.
            Parameters = obj.potentialFunc.Workspace.Parameters;
            fnames = fieldnames(Parameters);
            for i = 1:numel(fnames)
                value(i,1) = Parameters.(fnames{i})(1);
                stdev(i,1) = Parameters.(fnames{i})(2);
            end
            coeffTable = table(value,stdev,'RowNames',fnames,'VariableNames',{'Value','Std'});
        end
        
        function varargout = viewWorkspace(obj)
            % VIEWWORKSPACE Outputs the Workspace in a human-readable table.
            %
            % See also FUNC.VIEWWORKSPACE
            if nargout == 0
                obj.potentialFunc.viewWorkspace;
            else
                varargout{:} = obj.potentialFunc.viewWorkspace;
            end
        end
        
        function report(obj)
            % REPORT Displays a report of all the properties of the object
            %
            % See also ELECTROLYZERMODEL.VIEWWORKSPACE
            msg1 = ['\nElectrolyzer model properties:\n'...
                    '\n type: %s\n' ...
                    '\n electrolyte: %s\n' ...
                    '\n molarMassOfElectrolyte: %6.4f kg/mol\n'];
            fprintf(msg1,obj.type,obj.electrolyte,obj.molarMassOfElectrolyte)
            msg2 = ['\n potentialFunc:\n' ...
                    '\n   equation: %s\n' ...
                    '\n   Workspace:\n'];
            fprintf(msg2,obj.potentialFunc.equation)
            workspaceReport = obj.potentialFunc.viewWorkspace;
            if ~isempty(workspaceReport)
                disp(workspaceReport)
            else
                fprintf('      Empty Workspace\n')
            end
            fprintf('\n funcStorage:\n')
            if ~isempty(obj.funcStorage)
                disp(obj.funcStorage)
            else
                fprintf('   Empty funcStorage\n')
            end
            fprintf('\n PlottingCurves:')
            if ~isempty(obj.PlottingCurves)
                fprintf(' [%dx%d struct] \n',size(obj.PlottingCurves,1),size(obj.PlottingCurves,2))
            else
                fprintf(' No curves available. Perform model fit to get plotting curves.\n\n')
            end
        end
        
        
        function childObj = copy(obj)
            % COPY  Creates a full copy of the object with its own handle.
            %   Properties of the child object are no longer related to
            %   those of the parent.
            childObj = electrolyzerModel('type',obj.type,'electrolyte',obj.electrolyte);
            if ~func.isEmpty(obj.potentialFunc)
                childObj.addPotentials(obj.funcStorage.func.copy,'names',obj.funcStorage.name);
            end
            childObj.setInWorkspace(obj.potentialFunc.Workspace);
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
                                    'wtfrac',[])));
                    case 'NaOH'
                        obj.molarMassOfElectrolyte = 22.9898 + 15.9994 + 1.0079;
                        obj.setInWorkspace(struct(...
                                'Variables',struct(...
                                    'electrolyte',2,...
                                    'molarMassOfElectrolyte',obj.molarMassOfElectrolyte,...
                                    'molality',[],...
                                    'Molarity',[],...
                                    'wtfrac',[])));
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
        
        function synchronizeFuncStorage(obj)
            % SYNCHRONIZEFUNCSTORAGE  Synchronize Workspace structures in
            % the funcStorage.
            %  Modifies the Workspace structure of all the func objects in 
            %  the funcStorage to match the Workspace of the potentialFunc 
            %  object.
            for i = 1:height(obj.funcStorage)
                obj.funcStorage.func(i).replaceInWorkspace(obj.potentialFunc.Workspace,'rebuild')
            end
        end

        function recreatePotentialFunc(obj)
            % RECREATEPOTENTIALFUNC Recreate potentialFunc from the storage.
            obj.clearPotentials('rebuild');
            funcs = obj.funcStorage.func;
            obj.addPotentials(funcs,'rebuild');
        end
    end
    
end