classdef Experiment < handle
    %EXPERIMENT Defines a class of single experiment
    %   Onjects store data and perform simple functions for the data
    
    properties
        name='nimi';
        wid=10;
        volume;
        data
        timevector
        data_avg = timetable();
        fitresult
        UI
        
    end
    
    properties (SetAccess = private, GetAccess = private)
        len
        opts = delimitedTextImportOptions("NumVariables", 18);
        
   end
    
    methods
        function obj = Experiment(filename,start_time)
            %read data from file and create object
            obj.opts.DataLines = [24, Inf];
            obj.opts.Delimiter = "\t";
            obj.opts.VariableNames = ["time", "Voltage", "Current", "Flow_ele", "P1", "P2", "Ustack", "Udummy","Tcondenser", "Tin", "Tout", "Current_DMM6500", "Flush", "N2", "O2", "H2", "H2O", "BGA"];
            obj.opts.VariableTypes = ["double", "double", "double", "double", "double", "double", "double","double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];
            %obj.opts = setvaropts(obj.opts, 1, "InputFormat", "yyyy-MM-dd'T'HH:mm:ss");
            obj.opts.ExtraColumnsRule = "ignore";
            obj.opts.EmptyLineRule = "read";
            %obj.opts.VariableNamingRule = "preserve";
               
            obj.data = readtable(filename, obj.opts);
            obj.data.Time = datetime(start_time,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss") + seconds(obj.data.time);
            obj.data = table2timetable(obj.data,'Rowtimes','Time');
        end
        
        function calcSEC(obj,Nc)
            %calculate thoretical hydrogen production and SEC
            z = 2; %Number of electrons transferred per mole (-)
            F = 96485; %Faraday constant (As/mol) or (C/mol)
            M_H = 2.016;
            R_H = 89.9e-3; %Hydrogen density (kg/Nm^3)
            
           % Hydrogen production kg/s
           obj.data.H2prod_the = obj.data.Current/(z*F) * Nc * M_H / 1e3;
           obj.data.H2prod_meas = obj.data.Flush .* obj.data.H2/100 / 1000 / 60 * R_H;
    
           % SEC kWh/kg
           obj.data.SEC_the = obj.data.Current.*obj.data.Voltage/3.6e6 ./ obj.data.H2prod_the;
           obj.data.SEC_meas = obj.data.Current.*obj.data.Voltage/3.6e6 ./ obj.data.H2prod_meas;
    
        end
        
        function calcLoss(obj)
            cp_water = 4184; %Water specific heat capacity (J·kg−1⋅K−1), should be changed according to electrolyte
            obj.data.deltaT = obj.data.Tout-obj.data.Tin;
            obj.data.Ploss = cp_water .* obj.data.deltaT .* obj.data.Flow_ele/1e3/60;
        end
        
        function plotTemperature(obj,variablename)
           %plot temperature data
           figure()
            plot(obj.data.Time,eval(['obj.data.',variablename]));
            grid on
            xlabel('Time')
            ylabel('Temperature (\circC)')        
        end
         
         function plotVoltageCurrent(obj,name)
           %plot U and I as a function of time
           figure('name',name)
            yyaxis left
            plot(obj.data.Time,obj.data.Current);
            ylabel('Current(A)') 
            %ylim([0 12])
            grid on
            yyaxis right
            plot(obj.data.Time,obj.data.Voltage);
            %ylim([1 2.7])
            xlabel('Time')
            ylabel('Voltage (V)')        
         end
         
      function plotUIcurve(obj,name)
           %plot U and I as a function of time
           figure('name',name)
            %yyaxis left
            plot(obj.data_avg.Current,obj.data_avg.Voltage,'marker','x');
            xlabel('Current(A)') 
            grid on
            ylabel('Voltage (V)')        
         end
        
        function GetAveragedData(obj,timestrings,avg_time)
            %Take average data 
                %timestrings: end time
                %avg_time: averaged sample length (min)
            obj.timevector = datetime(timestrings,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
            obj.timevector(:,2) = obj.timevector - minutes(avg_time); %add start times for averaging
            
            for n=1:size(timestrings,1)
             tempdata = obj.data(timerange(obj.timevector(n,2),obj.timevector(n,1)),:); %create tempdata of the desired timerange
             obj.data_avg = [obj.data_avg;varfun(@mean,tempdata)]; %calculate the mean value of the tempdata
            end
            obj.data_avg.Time2 = obj.data_avg.Time+minutes(avg_time); %add new time vector with end times
            obj.data_avg.Properties.VariableNames=regexprep(obj.data_avg.Properties.VariableNames, 'mean_', ''); %remove "mean_" from variable name
        end
        
        function GetCurveFit(obj)
            xData = obj.data_avg.Current;
            yData = obj.data_avg.Voltage;

            xData(2) = [];
            yData(2) = [];
            
            % Set up fittype and options.
            ft = fittype( 'a + b*x + c*asinh(x/d)', 'independent', 'x', 'dependent', 'y' );
            fopts = fitoptions( 'Method', 'NoNlinearLeastSquares' );
            fopts.Display = 'Off';
            fopts.Lower = [0.3 0.0001 0 0.0]; %e0 k r
            fopts.StartPoint = [0.5 0.1 .1 .01];
            fopts.Upper = [.8 .1 1 .5];

            % Fit model to data.
            [obj.fitresult, gof] = fit( xData, yData, ft, fopts )
        end
        
        function GetUIdata(obj,Ivec)
           delta = 10e-3;
           ind(2) = 1;
            for n=1:length(Ivec)
                ind(1) = find(obj.data.Current(ind(2):end) > Ivec(n)-delta & obj.data.Current(ind(2):end) < Ivec(n)+delta,1) + ind(2)-1;
                ind(2) = min([find(obj.data.Current(ind(1):end) < Ivec(n)-delta | obj.data.Current(ind(1):end) > Ivec(n)+delta,1) , length(obj.data.Current(ind(1):end))])+ ind(1)-1;
                Time = obj.data.Time(ind(1):ind(2));
                Voltage = obj.data.Voltage(ind(1):ind(2));
                obj.UI.Time(1:length(Time),n) = Time;
                obj.UI.rTime(1:length(Time),n) = Time-Time(1);
                obj.UI.Voltage(1:length(Time),n) = Voltage;

            end
            obj.UI.Voltage(find(obj.UI.Voltage==0)) = NaN;        
         end
      
    end
end

