clear;
close all;
clc;


%% Creating function handle for fit

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pCat = 30; % bar
pAn = 2; % bar

load('JulichData.mat')

%% Create electrolyzer model objects
Emodel = electrolyzerModel('type','PEM');
Emodel.setParams(struct('Variables',struct('T',T,'pCat',pCat,'pAn',pAn)));
Emodel.addPotentials('ocv','act','ohm');

for j = 1
    
    [Imeas,sortInd] = sort(JulichUI(j).I);
    I = [Imeas JulichUI(j).Istd(sortInd)];
    U = [JulichUI(j).U(sortInd) JulichUI(j).Ustd(sortInd)];
    
    %% Fit
    
    % Particleswarm
    tic;
    [fitParams,gof] = Emodel.fitUI(U,I,'method','ps','weights','default');
    toc
    
    % Calculate voltage with calculate-method of func-object
    Ufit = Emodel.calculate('current',Imeas);
    
    % Plotting
    
    figure
    hold on;
    errorbar(I(:,1),U(:,1),U(:,2),U(:,2),I(:,2),I(:,2),'o')
    plot(Imeas,Ufit)
    xlabel("I (A)")
    ylabel("U (V)")
    legend("Data", "Fit", "Location", "Best")
    title(['Particleswarm: Jülich ' num2str(j)])
    hold off;
    
    %% Calculating temperature dependency of the different potentials
    
    T = (0:1:100)'+273.15;
    pCat = 1.15; % bar
    pAn = 1.15; % bar
    
    Emodel.replaceParams('T',T,'pCat',pCat,'pAn',pAn,'current',2.2);
    
    Urev = reversible(T,6).calculate;
    Uocv = Emodel.funcStorage.func(1).calculate;
    Uact = Emodel.funcStorage.func(2).calculate;
    Uohm = Emodel.funcStorage.func(3).calculate;
    
    figure
    plot(T,[Urev Uocv Uocv+Uact Uocv+Uact+Uohm])
    xlabel('Temperature')
    ylabel('Voltage')
    legend('Urev','Uocv','Uocv + Uact','Uocv + Uact + Uohm')
end