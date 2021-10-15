clear;close all;clc

%%
ePEM = electrolyzerModel('type','pem');

T = 273.15 + (30:10:60)';

Workspace.Variables.T = T;
Workspace.Variables.pCat = 20;
Workspace.Variables.pAn = 2;

Workspace.Coefficients.alpha = 0.5;
Workspace.Coefficients.j0 = 1e-12;
Workspace.Coefficients.r = 0.1;

ePEM.setParams(Workspace)

ePEM.addPotentials('nernst','ohmic','activation')

current = (0.02:0.001:1.5);

voltage = ePEM.calculate('current',current);

figure
plot(current,voltage)
legend('30','40','50','60')
 
%% Separate potentials
Uocv = nernst('pem');
Uocv.replaceParams(Workspace)

figure
plot(current,ones(size(current)).*Uocv.calculate())
title('Uocv')
xlabel('Current density (A/cm^2)')
ylabel('Voltage (V)')
legend('30','40','50','60')

Urev = nernstReversible(1);
Urev.replaceParams(Workspace)

figure
plot(current,ones(size(current)).*Urev.calculate())
title('Urev')
xlabel('Current density (A/cm^2)')
ylabel('Voltage (V)')
legend('30','40','50','60')

Upres = nernstPressureCorrection('pem');
Upres.replaceParams(Workspace)

figure
plot(current,ones(size(current)).*Upres.calculate())
title('Upres')
xlabel('Current density (A/cm^2)')
ylabel('Voltage (V)')
legend('30','40','50','60')

Uact = activation();
Uact.replaceParams(Workspace)

figure
plot(current,Uact.calculate('current',current)')
title('Uact')
xlabel('Current density (A/cm^2)')
ylabel('Voltage (V)')
legend('30','40','50','60')

Uohm = ohmic();
Uohm.replaceParams(Workspace)

figure
plot(current,Uohm.calculate('current',current)')
title('Uohm')
xlabel('Current density (A/cm^2)')
ylabel('Voltage (V)')


%% Test Butler-Volmer
[F,R,n_e] = getConstants();
j0ref = 1e-8;
Tref = T(1);
Ea = (76+18)*1e3; % J/mol, activation energy
gammaM = 100; % Roughness factor
a = 0.5;
UactperTf = (-15:0.01:15);
Uact = UactperTf*R.*T/(F*n_e);
j0 = gammaM*exp(Ea/R*(1/Tref-1./T))*j0ref;
j = j0.*(exp(a.*UactperTf) - exp((a-1).*UactperTf));

figure
plot((j.*(abs(j)<1.5))',(Uact.*(abs(j)<1.5))')
title('ButtlerVolmer')
xlabel('Current density (A/cm^2)')
ylabel('Voltage (V)')
legend('30','40','50','60','location','best')










