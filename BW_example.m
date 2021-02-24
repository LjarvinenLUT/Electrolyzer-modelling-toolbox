clear;
close all;
addpath("Utils");

% Script for plotting different models for activity overpotential
[fileOutputPath,FigSizeWidth,FigSizeHeight] = setFigureParameters();
%% Global parameters
N_A = 6.02214076e23; % 1/mol, Avogadro constant
k_B = 1.380649e-23; % J/K, Boltzmann constant
Q_e = 1.602176634e-19; % C, Elemental charge
global F n_e R
F = Q_e*N_A; % C/mol, Faraday constant
R = k_B*N_A; % J/K/mol, Universal gas constant
n_e = 2; % Number of electrons transferred in one reaction

f = R/(n_e*F);

%% Drawing of Buttler-Volmer and other activity over potential approximations

%% Buttler-Volmer with variable alpha

UactperfT = -10:0.01:10;
alpha = (0.3:0.1:0.7)';
jperj0 = exp(alpha.*UactperfT)-exp((alpha-1).*UactperfT);


fig = figure('Name','Buttler-Volmer for variable alpha');
hold on;
for i = 1:numel(alpha)
    p(i) = plot(jperj0(i,:),UactperfT);
end
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j/j_0$')
xlim([-20 20])
legendflex([p(1) p(2) p(3) p(4) p(5)],{'0.3','$0.4$','$0.5$','$0.6$','0.7'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10],'title','$\alpha$')

fullFileOutput = fileOutputPath+"BV_var_alpha.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');

%% Buttler-Volmer compared to Tafel

jperj0Tafel1 = exp(alpha.*UactperfT);
jperj0Tafel2 = -exp((alpha-1).*UactperfT);

fig = figure('Name','Buttler-Volmer vs. Tafel for alpha=1/2');
hold on;
p1 = plot(jperj0(3,:),UactperfT);
p2 = plot(jperj0Tafel1(3,:),UactperfT);
p3 = plot(jperj0Tafel2(3,:),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j/j_0$')
xlim([-20 20])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BV_vs_Tafel.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');

fig = figure('Name','Logarithmic Buttler-Volmer vs. Tafel for alpha=1/2');
hold on;
p1 = plot(log(abs(jperj0(3,:))),UactperfT);
p2 = plot(log(abs(jperj0Tafel1(3,:))),UactperfT);

p3 = plot(log(abs(jperj0Tafel2(3,:))),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$\ln(|j/j_0|)$')
% xlim([-60 60])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BV_vs_Tafel_log.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');

<<<<<<< .merge_file_a12580

fig = figure('Name','Logarithmic Buttler-Volmer vs. Tafel for alpha=0.7');
=======
f = figure('Name','Logarithmic Buttler-Volmer vs. Tafel for alpha=0.7');
>>>>>>> .merge_file_a09848
hold on;
p1 = plot(log(abs(jperj0(5,:))),UactperfT);
p2 = plot(log(abs(jperj0Tafel1(5,:))),UactperfT);
p3 = plot(log(abs(jperj0Tafel2(5,:))),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$\ln(|j/j_0|)$')
% xlim([-60 60])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BV_vs_Tafel_log_alpha_0.7.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');

%% Buttler-Volmer compared to asinh approximation

jperj0sinh = 2*sinh(alpha.*UactperfT);

fig = figure('Name','Buttler-Volmer vs sinh approximation with variable alpha');
hold on;
colorOrder = get(gca, 'ColorOrder');
p1 = [];
p2 = [];
for i = 1:2:numel(alpha)
    j = i-(i-1)/2;
    p1(j) = plot(jperj0(i,:),UactperfT,'--','Color',colorOrder(j,:));
    p2(j) = plot(jperj0sinh(i,:),UactperfT,'Color',colorOrder(j,:));
end
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j/j_0$')
xlim([-10 10])
legendflex([p2(1) p2(2) p2(3)],{'0.3','0.5','0.7'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10],'title','$\alpha$')

fullFileOutput = fileOutputPath+"BV_vs_sinh_var_alpha.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');


%% asinh approximation for combined electrode

i = 3;
alphac = 1
alphaa = 1
j0c = 1e-3
j0a = j0c*2
j = jperj0(i,:)*j0c;
jperj0c = j/j0c;
jperj0a = j/j0a;

UactperfTc = 1/alphac*asinh(jperj0c/2);
UactperfTa = 1/alphaa*asinh(jperj0a/2);

UactperfTtot = UactperfTc + UactperfTa;

% Fit combined sinh function
fitfun = fittype(@(fitalpha,fitj0,x) 1/fitalpha*asinh(x/fitj0));
[fitted_curve,gof] = fit(j',UactperfTtot',fitfun,'StartPoint',[min(alphac,alphaa) min(j0c,j0a)]);

coeffvals = coeffvalues(fitted_curve);

alphatot = coeffvals(1)
j0tot = coeffvals(2)

figure('Name','combined sinh approximation')
plot(j,UactperfTc,j,UactperfTa,j,UactperfTtot,j,fitted_curve(j),'--')
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j$')
xlim([0 max(j)])

%% Three approximations compared

fig = figure('Name','Buttler-Volmer vs sinh approximation vs Tafel with alpha = 0.3');
hold on;
p1 = plot(jperj0(1,:),UactperfT);
p2 = plot(jperj0sinh(1,:),UactperfT);
p3 = plot(jperj0Tafel1(1,:),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j/j_0$')
xlim([0 10])
legendflex([p1 p2 p3],{'Buttler-Volmer','sinh approximation','Tafel approximation'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10],'title','$\alpha = 0.3$')

fullFileOutput = fileOutputPath+"BV_vs_sinh_vs_tafel.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');
