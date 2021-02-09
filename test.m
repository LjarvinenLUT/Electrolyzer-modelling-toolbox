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

T = linspace(273,373,100)';
U = [];
figure
hold on;
for i = 1:6
    U_0 = reversible(T,i);
    U = [U U_0];
    plot(T,U_0)
end
legend('1','2','3','4','5','6')
hold off;

difU = abs(U-mean(U,2));

figure
plot(difU)