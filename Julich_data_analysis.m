%% Script for Jülich data analysis
close all; clear; clc;

filename = {"Julich_Uj_data1.xlsx","Julich_Uj_data2.xlsx","Julich_Uj_data3.xlsx"};

JulichUI = struct('U',[],'Ustd',[],'I',[],'Istd',[]);


for j = 1:3
    
    %% Read data from file
    opts = detectImportOptions(filename{j});
    opts.SelectedVariableNames = 1:3;
    if j == 3
        opts.DataRange = '1:34990';
    end
    M = readmatrix(filename{j},opts);
    
    %% Timeseries, removing NaNs
    nans = any(isnan(M),2); % Find and remove any nan values from the vectors
    
    t = M(~nans,1);
    Umeas = M(~nans,2);
    Imeas = M(~nans,3);
    [Imeassort,I] = sort(Imeas);
    Umeassort = Umeas(I);
    
    figure
    yyaxis left
    scatter(1:length(Imeas),Imeas)
    ylabel('I')
    yyaxis right
    scatter(1:length(Umeas),Umeas)
    ylabel('U')
    title('Measured Jülich UI timeseries')
    xlabel('timestep')
    
    % Moving average filtering for j == 3
    if j == 3
        M = 50;
        h = 1/M*ones(1,M);
        Ifiltered = filter(h,1,Imeas);
        
        figure
        yyaxis left
        scatter(1:length(Ifiltered),Ifiltered)
        ylabel('I')
        yyaxis right
        scatter(1:length(Umeas),Umeas)
        ylabel('U')
        title('Filtered Jülich UI timeseries')
        xlabel('timestep')
    else
        Ifiltered = Imeas;
    end
    
    
    

    
    %% Difference between adjacent measurements
    difU = diff(Umeas);
    difI = diff(Ifiltered);
    
    figure
    yyaxis left
    scatter(1:length(difI),difI)
    ylabel('dI')
    yyaxis right
    scatter(1:length(difU),difU)
    ylabel('dU')
    title('Jülich UI timeseries difference')
    xlabel('timestep')
    
    if j == 1
        difthresh = 1e-3;
    elseif j == 2
        difthresh = 6e-3;
    elseif j == 3
        difthresh = [1e-2*ones(18300,1);0.03*ones(length(difI)-18300,1)];
    end
    Ilevel = [abs(difI)<difthresh;false];
    
    figure
    yyaxis left
    scatter(1:length(Imeas),Imeas)
    ylabel('I')
    yyaxis right
    plot(Ilevel)
    ylabel('dI<threshold')
    title('Jülich UI timeseries difference')
    xlabel('timestep')
    
    figure
    yyaxis left
    scatter(1:length(Umeas),Umeas)
    ylabel('U')
    yyaxis right
    plot(Ilevel)
    ylabel('dI<1e-4')
    title('Jülich UI timeseries difference')
    xlabel('timestep')
    
    
    %% Pick UI values
    
    prev_level = false;
    I = [];
    U = [];
    Istd = [];
    Ustd = [];
    Itemp = [];
    Utemp = [];
    
    for i = 1:length(t)
        if ~Ilevel(i)
            if prev_level % If level current ends
                I = [I;mean(Itemp)];
                Istd = [Istd;sqrt(var(Itemp))];
                difUtemp = abs(diff(Utemp));
                Utemplevel = difUtemp<3e-4;
                U = [U;mean(Utemp(Utemplevel))];
                Ustd = [Ustd;sqrt(var(Utemp(Utemplevel)))];
                Itemp = [];
                Utemp = [];
            end
            prev_level = false;
        elseif Ilevel(i)
            Itemp = [Itemp Imeas(i)];
            Utemp = [Utemp Umeas(i)];
            prev_level = true;
        end
    end
    
    nans = isnan(U);
    I = I(~nans);
    Istd = Istd(~nans);
    U = U(~nans);
    Ustd = Ustd(~nans);
    JulichUI(j).U = U;
    JulichUI(j).I = I;
    JulichUI(j).Ustd = Ustd;
    JulichUI(j).Istd = Istd;
    
    %% Final measured UI curve
    figure
    errorbar(I,U,Ustd,Ustd,Istd,Istd,'o')
    title('Jülich UI measurement')
    xlabel('I')
    ylabel('U')

end

%% Saving

save('JulichData.mat','JulichUI')