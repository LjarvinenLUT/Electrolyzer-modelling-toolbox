clc;
clear all;
close all;

%% Set figure properties
SetFigureProperties()

%% Test time string write
%  apu = datetime('2021-09-07T08:21:00','InputFormat',"yyyy-MM-dd'T'HH:mm:ss");
%  apu2 = apu + minutes([0:17])
%  timestrings = datestr(apu2,'yyyy-mm-ddTHH:MM:ss')
% timestrings2 = ['2021-09-08T08:21:00';'2021-09-08T08:22:00';'2021-09-08T08:23:00';'2021-09-08T08:24:00';'2021-09-08T08:25:00';'2021-09-08T08:26:00';'2021-09-08T08:27:00';'2021-09-08T08:28:00';'2021-09-08T08:29:00';...
%                 '2021-09-08T08:30:00';'2021-09-08T08:31:00';'2021-09-08T08:32:00';'2021-09-08T08:33:00';'2021-09-08T08:34:00';'2021-09-08T08:35:00';'2021-09-08T08:36:00';'2021-09-08T08:37:00';'2021-09-08T08:38:00'];
% jep = datetime(timestrings,'InputFormat',"yyyy-MM-dd'T'HH:mm:ss");            
%  a
%% Create objects for each experiment
for n=1:4
    
    switch n
         case 1
            filename = 'IUCURVE_30C.lvm';
            start_time = '2021-09-07T13:59:22';
            timestrings = ['2021-09-07T14:03:00';'2021-09-07T14:06:30';'2021-09-07T14:28:30';'2021-09-07T14:44:00';'2021-09-07T15:00:00'...
                ;'2021-09-07T15:16:00';'2021-09-07T15:31:40';'2021-09-07T15:32:30';'2021-09-07T15:33:45';'2021-09-07T15:34:40'...
                ;'2021-09-07T15:35:30';'2021-09-07T15:36:10';'2021-09-07T15:36:50';'2021-09-07T15:37:30';'2021-09-07T15:38:10';'2021-09-07T15:38:50'];
        case 2
            filename = 'Slow_08-09-2021_60C.lvm';
            start_time = '2021-09-08T08:20:56';
            apu = datetime('2021-09-08T08:21:00','InputFormat',"yyyy-MM-dd'T'HH:mm:ss"); apu2 = apu + minutes([0:17]); timestrings = datestr(apu2,'yyyy-mm-ddTHH:MM:ss')
            %timestrings = ['2021-09-08T08:21:00';'2021-09-08T08:22:00';'2021-09-08T08:23:00';'2021-09-08T08:24:00';'2021-09-08T08:25:00';'2021-09-08T08:26:00';'2021-09-08T08:27:00';'2021-09-08T08:28:00';'2021-09-08T08:29:00';...
              %  '2021-09-08T08:30:00';'2021-09-08T08:31:00';'2021-09-08T08:32:00';'2021-09-08T08:33:00';'2021-09-08T08:34:00';'2021-09-08T08:35:00';'2021-09-08T08:36:00';'2021-09-08T08:37:00';'2021-09-08T08:38:00'];
            
        case 3 %NO BGA results
            filename = 'Slow_06-09-2021_30C.lvm';
            start_time = '2021-09-06T11:37:26';
            timestrings = ['2021-09-06T11:40:00';'2021-09-06T11:57:00';'2021-09-06T12:20:00'];
        case 4 %NO BGA results
            filename = 'Slow_06-09-2021_60C.lvm';
            start_time = '2021-09-06T12:28:07';
            timestrings = ['2021-09-06T12:39:00';'2021-09-06T12:59:00';'2021-09-06T13:15:00';'2021-09-06T13:33:00';'2021-09-06T13:49:00'];
    end
    
    avg_time = 0.25; %averaging time (min)
    Nc = 10; %number of cells in series
    
    % make object with data from file
    e(n) = Experiment(filename,start_time)
    
    calcSEC(e(n),Nc)
    calcLoss(e(n))
    % plot desired temperature
    plotTemperature(e(n),'Tin')
    
    plotVoltageCurrent(e(n),'Voltage & Current')
    
    % get averaged data at selected time instants
    GetAveragedData(e(n),timestrings,avg_time)
    
    plotUIcurve(e(n),'UI')
    
end

%%

n=2;
figure('name','Time series')
    yyaxis left
    plot(e(n).data.Time,e(n).data.Current);
    hold on
    plot(e(n).data_avg.Time,e(n).data_avg.Current,'marker','x','linestyle','none');
    ylabel('Current(A)')
    %ylim([0 12])
    grid on
    yyaxis right
    plot(e(n).data.Time,e(n).data.Voltage);
        hold on
    plot(e(n).data_avg.Time,e(n).data_avg.Voltage,'marker','x','linestyle','none');
    %ylim([1 2.7])
    xlabel('Time')
    ylabel('Voltage (V)')

figure('name','UI-comparison')
for n=1:4
    plot(e(n).data_avg.Current,e(n).data_avg.Voltage,'marker','x');
    hold on
end
xlabel('Current(A)')
grid on
ylabel('Voltage (V)')
legend('30 °C (8.9.2021)','60 °C (8.9.2021)','30 °C (6.9.2021)','60 °C (6.9.2021)','location','southeast')

figure('name','Loss power')
for n=1:4
    plot(e(n).data_avg.Current,e(n).data_avg.Ploss,'marker','x');
    hold on
end
xlabel('Current(A)')
grid on
ylabel('Loss power (W)')
legend('30 °C','60 °C','location','southeast')

%%
n=2
figure('name','H2 content')
plot(e(n).data_avg.Current,e(n).data_avg.H2,'marker','x')
hold on
plot(e(n).data_avg.Current,e(n).data_avg.BGA,'marker','x')
xlabel('Current(A)')
grid on
ylabel('Hydrogen content (%)')
legend('MS','BGA')

%%
n=2
figure('name','Flow')
plot(e(n).data.Time,e(n).data.Flush,'marker','x')
xlabel('Time (-)')
grid on
ylabel('Air flow (L/min)')

%% Test modelling toolbox
eModel = electrolyzerModel('type','alkaline','electrolyte','NaOH');
T = 273.15 + 60; % Temperature in kelvin
ps = 1; % System pressure in bara
wtfrac = 20; % concentration as weight percentage
m = eModel.wtfrac2mol(wtfrac)
Workspace = struct('Variables',struct('T',T,'ps',ps,'m',m))
func.isWorkspace(Workspace)
eModel.setParams(Workspace)
eModel.viewWorkspace;
eModel.addPotentials('nernst','ohmic','activation')

method = "PS";
weights = "none";
voltageData = e(1).data_avg.Voltage/10;
currentData = e(1).data_avg.Current/(pi*2.5^2);
[fitParams,gof] = eModel.fitUI(voltageData,currentData,'method',method,'weights',weights);
eModel.showUI
disp(fitParams)