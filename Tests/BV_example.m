clear;
close all;

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
xlim([-5 5])
ylim([-10 10])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BV_vs_Tafel_log.pdf";
% saveas(gcf,fullFileOutput);
% export_fig(gcf, fullFileOutput);
fig.PaperUnits = 'centimeters';
fig.PaperSize = [FigSizeWidth FigSizeHeight]; % Make the "page" just big enough to hold the size output I want
% f.PaperPosition(1:2) = [0 0]; % Specify that we start at the lower-left corner of the page
fig.PaperPositionMode = 'auto';
print(fullFileOutput,'-dpdf','-r0');


fig = figure('Name','Logarithmic Buttler-Volmer vs. Tafel for alpha=0.7');
hold on;
p1 = plot(log(abs(jperj0(5,:))),UactperfT);
p2 = plot(log(abs(jperj0Tafel1(5,:))),UactperfT);
p3 = plot(log(abs(jperj0Tafel2(5,:))),UactperfT);
hold off;
yline(0)
ylabel('$U_{\mathrm{act}}/fT$')
xline(0)
xlabel('$\ln(|j/j_0|)$')
xlim([-5 5])
ylim([-10 10])
legendflex([p1 p2 p3],{'Buttler-Volmer','Tafel equation for positive current','Tafel equation for negative current'},'xscale', 0.4,'anchor',{'se','se'},'buffer',[-10 10])

fullFileOutput = fileOutputPath+"BV_vs_Tafel_log_alpha_0_7.pdf";
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
alphac = 0.1:0.1:0.9;
alphaa = alphac;
j0c = linspace(0.1,10,50)*1e-3;
j0a = j0c;
j = (jperj0(i,:))'*1e-3;
jperj0c = j./j0c;
jperj0a = j./j0a;

filename = 'alpha_j0_simulation_results.mat';
if exist(filename,'file') == 2
    load(filename)
else

    fitfun = fittype(@(fitalpha,fitj0,x) 1/fitalpha*asinh(x/fitj0));
    
    for i = 1:length(alphac)
        for k = 1:length(alphaa)
            for t = 1:length(j0c)
                for s = 1:length(j0a)
                    
                    UactperfTc = 1/alphac(i)*asinh(jperj0c(:,t)/2);
                    UactperfTa = 1/alphaa(k)*asinh(jperj0a(:,s)/2);
                    
                    UactperfTtot = UactperfTc + UactperfTa;
                    
                    % Fit combined sinh function
                    [fitted_curve,gof] = fit(j,UactperfTtot,fitfun,'StartPoint',[min(alphac(i),alphaa(k)) min(j0c(t),j0a(s))]);
                    
                    coeffvals = coeffvalues(fitted_curve);
                    
                    alphatot(i,k,t,s) = coeffvals(1);
                    j0tot(i,k,t,s) = coeffvals(2);
                    
                end
            end
        end
    end
    
    save(filename,'alphatot','j0tot')
end

%% alphatot as a function of alpha1 and alpha2
figure('Name','alpha as a function of alpha_a and alpha_c, j_0a and j_0c = 4.9e-3')
surf(alphac,alphaa,alphatot(:,:,25,25))

figure('Name','alpha as a function of alpha_a and alpha_c, j_0a and j_0c = 0.1e-3')
surf(alphac,alphaa,alphatot(:,:,1,1))

alpha_funkalpha = (alphac.*alphaa')./(alphac + alphaa');

figure('Name','alpha1*alpha2/(alpha1 + alpha2)')
surf(alphac,alphaa,alpha_funkalpha)

%% alphatot as a function of j01 and j02
figure('Name','alpha as a function of j_0a and j_0c, alpha_a and alpha_c = 1/2')
surf(j0c,j0a,(permute(alphatot(5,5,:,:),[3 4 1 2])/alpha_funkalpha(5,5)))

figure('Name','alpha as a function of j_0a and j_0c, alpha_a and alpha_c = 0.9')
surf(j0c,j0a,(permute(alphatot(9,1,:,:),[3 4 1 2])/alpha_funkalpha(9,1)))

%% j0tot as a function of alpha1 and alpha2
figure('Name','j_0 as a function of alpha_a and alpha_c, j_0a and j_0c = 4.9e-3')
surf(alphac,alphaa,j0tot(:,:,25,25))

%% j0tot as a function of j01 and j02
figure('Name','j_0 as a function of j_0a and j_0c, alpha_a and alpha_c = 1/2')
surf(j0c,j0a,permute(j0tot(5,5,:,:),[3 4 1 2]))


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
