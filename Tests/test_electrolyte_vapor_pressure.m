% KOH electrolyte water vapor pressure comparison

[fileOutputPath,FigSizeWidth,FigSizeHeight] = setFigureParameters();

m = 0:0.1:18; % Molality, mol/kg of solvent
T = 313.15; % Temperature

MH2O = (2*1.008 + 16)*1e-3; % kg/mol, molar mass of water

molfracSolute = m*MH2O; % Molar fraction of solute
molfracH2O = 1-molfracSolute; % Molar fraction of H2O


[psvEl1,~] = electrolyte_parameters(T,m,'KOH',1); % Experimental fit
[psvEl2,~] = electrolyte_parameters(T,m,'KOH',2); % Ideal solution approximation
psvElMeas = 10.^[-1.13224 -1.16941 -1.22185 -1.29073 -1.37469 -1.47108 -1.57675 -1.69250 -1.81816 -1.94310];
mMeas = [0 2 4 6 8 10 12 14 16 18];

figure
plot(molfracH2O,psvEl1,molfracH2O,psvEl2)
ylabel('Saturated vapor pressure of water, bar','interpreter','latex')
xlabel('Molar fraction of water in the solution','interpreter','latex')
legend('Experimental fit','Ideal solution approximation')

fig = figure;
plot(m,psvEl1,m,psvEl2)
hold on
scatter(mMeas,psvElMeas)
ylabel('Saturated vapor pressure of water, bar','interpreter','latex')
xlabel('KOH molality, mol/kg of solvent','interpreter','latex')
legend('Experimental fit','Ideal solution approximation','Measured')
title('Water vapor pressure in KOH solution, T = 313.15 K')

fullFileOutput = fileOutputPath+"saturated_vapor_pressure_KOH_solution.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');
