close all;
clear;
clc;

%% Global parameters
N_A = 6.02214076e23; % 1/mol, Avogadro constant
k_B = 1.380649e-23; % J/K, Boltzmann constant
Q_e = 1.602176634e-19; % C, Elemental charge

F = Q_e*N_A; % C/mol, Faraday constant
R = k_B*N_A; % J/K/mol, Universal gas constant
n_e = 2; % Number of electrons transferred in one reaction

f = R/(n_e*F);

Umax = 1;
fT = f*(273.15+25);

UactperfT = linspace(-Umax/fT,Umax/fT,100);
alpha = (0.1:0.1:0.9)';
jperj0 = exp(alpha.*UactperfT)-exp((alpha-1).*UactperfT);

j0(1,1,:) = logspace(-9,-3,7);

j = jperj0.*j0;

figure
plot(UactperfT,permute(j(5,:,:),[3 2 1]))
xlim([-40 40])
ylim([-1 1])
ylabel('Current density')
xlabel('$U/fT$')

%%

x = 0:0.01:10;
a = logspace(-1,2,7)';
y = a.*exp(x);

figure
plot(x,y)
xlim([0 6])
ylim([0 300])
ylabel('y')
xlabel('x')

y2 = linspace(min(y,[],'all'),300,100);
x2 = log(y2./a);

figure
plot(y2,x2)
ylim([0 6])
xlim([0 300])
ylabel('x')
xlabel('y')