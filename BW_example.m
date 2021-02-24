clear;
close all;

% Script for plotting different models for activity overpotential
[fileOutputPath,FigSizeWidth,FigSizeHeight] = setFigureParameters;

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


f = figure('Name','Buttler-Volmer for variable alpha');
hold on;
for i = 1:numel(alpha)
    p(i) = plot(jperj0(i,:),UactperfT);
end
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j/j_0$')
xlim([-60 60])
legendflex([p(1) p(2) p(3) p(4) p(5)],{'0.3','$0.4$','$0.5$','$0.6$','0.7'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10],'title','$\alpha$')

fullFileOutput = fileOutputPath+"BW_variable_alpha.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
f.PaperUnits = 'centimeters';
f.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
f.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');

%% Buttler-Volmer compared to other approximations

jperj0Tafel1 = exp(alpha.*UactperfT);
jperj0Tafel2 = -exp((alpha-1).*UactperfT);

f = figure('Name','Buttler-Volmer vs. Tafel for alpha=1/2');
hold on;
p1 = plot(jperj0(3,:),UactperfT);
p2 = plot(jperj0Tafel1(3,:),UactperfT);
p3 = plot(jperj0Tafel2(3,:),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$j/j_0$')
xlim([-30 30])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BW_vs_Tafel.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
f.PaperUnits = 'centimeters';
f.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
f.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');

f = figure('Name','Logarithmic Buttler-Volmer vs. Tafel for alpha=1/2');
hold on;
p1 = plot(log(abs(jperj0(3,:))),UactperfT);
p2 = plot(log(abs(jperj0Tafel1(3,:))),UactperfT);

p3 = plot(log(abs(jperj0Tafel2(3,:))),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$\ln(j/j_0)$')
% xlim([-60 60])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BW_vs_Tafel_log.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
f.PaperUnits = 'centimeters';
f.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
f.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');


f = figure('Name','Logarithmic Buttler-Volmer vs. Tafel for alpha=0.7');
hold on;
p1 = plot(log(abs(jperj0(5,:))),UactperfT);
p2 = plot(log(abs(jperj0Tafel1(5,:))),UactperfT);

p3 = plot(log(abs(jperj0Tafel2(5,:))),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$\ln(j/j_0)$')
% xlim([-60 60])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BW_vs_Tafel_log_alpha_0.7.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
f.PaperUnits = 'centimeters';
f.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
f.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');