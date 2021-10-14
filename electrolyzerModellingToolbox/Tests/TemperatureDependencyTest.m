ePEM = electrolyzerModel('type','pem');

Workspace.Variables.T = 273.15 + (30:10:60)';
Workspace.Variables.pCat = 20;
Workspace.Variables.pAn = 2;

Workspace.Coefficients.alpha = 0.3;
Workspace.Coefficients.j0 = 1e-4;
Workspace.Coefficients.r = 0.1;

ePEM.setParams(Workspace)

ePEM.addPotentials('nernst','ohmic','activation')

current = (0.02:0.001:1.5);

voltage = ePEM.calculate('current',current);

figure
plot(current,voltage)
legend('30','40','50','60')







