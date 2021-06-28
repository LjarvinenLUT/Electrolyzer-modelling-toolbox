
% Data import
filename = 'alkali_data_export.csv';

opts = detectImportOptions(filename);
opts.VariableTypes{1} = 'char';
opts.SelectedVariableNames = 1;
dt = datetime(readmatrix(filename,opts));
opts.SelectedVariableNames = 2:3;
M = readmatrix(filename,opts);
nans = any(isnan(M),2);
U = M(~nans,1);
I = M(~nans,2);
dt = dt(~nans);

data = {dt,U,I};

%% Pressures import
filename = 'alkali_pressure_export.csv';

opts = detectImportOptions(filename);
opts.VariableTypes{1} = 'char';
opts.SelectedVariableNames = 1;
dt = datetime(readmatrix(filename,opts));
opts.SelectedVariableNames = 2:3;
P = readmatrix(filename,opts);
nans = any(isnan(P),2);
P = mean(P(~nans),2);
dt = dt(~nans);

Pdata = {dt,P};

%% Temperature import
filename = 'alkali_temperature_export.csv';

opts = detectImportOptions(filename);
opts.VariableTypes{1} = 'char';
opts.SelectedVariableNames = 1;
dt = datetime(readmatrix(filename,opts));
opts.SelectedVariableNames = 2:4;
T = readmatrix(filename,opts);
nans = any(isnan(T),2);
T = mean(T(~nans),2);
dt = dt(~nans);

Tdata = {dt,T};


%% Extracting average information from correct times
times = {'10:27','10:37';...
    '10:52','11:02';...
    '11:17','11:27';...
    '11:55','12:05';...
    '12:28','12:38';...
    '13:26','13:36';...
    '14:07','14:17';...
    '14:39','14:49';...
    '15:10','15:20';...
    '15:28','15:38';...
    '16:05','16:15'}; % Measurement times
measurement_times = datetime(strcat('2020-04-06',{' '},times));

U = nan(length(measurement_times(:,1)),2);
I = nan(size(U));
T = nan(size(U));
P = nan(size(U));

for i = 1:length(measurement_times(:,1))
    % Voltage and current
    t = data{:,1}>=measurement_times(i,1)&data{:,1}<=measurement_times(i,2);
    temp = cell2mat(data(:,2:3));
    U(i,:) = [mean(temp(t,1)) std(temp(t,1))];
    I(i,:) = [mean(temp(t,2)) std(temp(t,2))];
    % Pressure
    t = Pdata{:,1}>=measurement_times(i,1)&Pdata{:,1}<=measurement_times(i,2);
    temp = Pdata{:,2};
    P(i,:) = [1+mean(temp(t)) std(temp(t))]; % Conversion from barg to bara
    % Temperature
    t = Tdata{:,1}>=measurement_times(i,1)&Tdata{:,1}<=measurement_times(i,2);
    temp = Tdata{:,2};
    T(i,:) = [273.15+mean(temp(t)) std(temp(t))]; % Conversion from degC to K
end

A = pi*5^2;        % effective area of a single cell [cm^2]
N = 30;         % number of cells

j = I./A;
U = U./N;

AlkaliUI = struct('U',U,'j',j,'T',T,'P',P);

save('AlkaliData.mat','AlkaliUI')


