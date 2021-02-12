% 
close all;
clear;

N_A = 6.02214076e23; % 1/mol, Avogadro constant
k_B = 1.380649e-23; % J/K, Boltzmann constant
Q_e = 1.602176634e-19; % C, Elemental charge
global F n_e R
F = Q_e*N_A; % C/mol, Faraday constant
R = k_B*N_A; % J/K/mol, Universal gas constant
n_e = 2; % Number of electrons transferred in one reaction

% Measured, from CRC Handbook of Chemistry and Physics
Tmeas = [298.15;300;373.21];
deltaG = [237.141;236.839;225.160]*1e3;
Umeas = deltaG/(n_e*F);

% Modelled
T = sortrows([linspace(273,373,101) 298.15 373.21]');
U = [];
difU = [];
figure
hold on;
for i = 1:6
    U_0 = reversible(T,i);
    U = [U U_0];
    plot(T,U_0)
end
scatter(Tmeas,Umeas,'x')
hold off;
legend('1','2','3','4','5','6','tabled')
title('Reversible voltage')
xlabel('Temperature [K]')
ylabel('Voltage [V]')


% Difference of measured to modelled
ComInd = [find(T==Tmeas(1));find(T==Tmeas(2));find(T==Tmeas(3))];
UforComp = U(ComInd,:);
difU = (UforComp-Umeas);


figure
hold on
for i = 1:length(difU(1,:))
    scatter(Tmeas,difU(:,i))
end
yline(0,'--')
hold off
legend('1','2','3','4','5','6')
title('Difference to measured value')
xlabel('Temperature [K]')
ylabel('Difference in voltage [V]')