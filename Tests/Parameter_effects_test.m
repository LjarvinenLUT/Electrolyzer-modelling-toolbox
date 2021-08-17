%% Script for visualising effects of different parameters in the UI curve fit
close all; clear; clc;

%% Global parameters
[F,R,n_e] = getConstants();

f = R/(n_e*F);


%% Creating function handle for fit

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pH2 = 30; % bar
pO2 = 2; % bar
Uocv = nernst(T,pH2,pO2,'type',type);
Uact = activation('model',2);
Uohm = ohmic();
Ucon = concentration();

Uforfit = @(j0,alpha,r,jL,Uerr,Current) Uocv + Uact(j0,alpha,Current) + Uohm(r,Current) + Ucon(jL,Current) + Uerr;

%% Creating test data

% Starting values for parameters
alpha = [0.25 0.3 0.75]; %
j0 = [1.5e-5 1.5e-6 1.5e-7]; % A/cm^2, exchange current density
r = [0.5 1 1.5]; % Ohm, total resistance
jL = [1 1.5 2]; % A/cm^2, limiting current density
Uerr = [-0.1 0 0.1]; % V, constant voltage error

parameters = [alpha;j0;r;jL;Uerr];

parameterNames = {'alpha','j0','r','jL','Uerr'};

U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits

for i = 1:length(parameters)

    Umeas = [[] [] []]; % "Measured" voltage
    jmeas = [[] [] []]; % "Measured" current based on Buttler-Volmer equation

    for j = 1:3

        useParameters = parameters(:,2);
        useParameters(i) = parameters(i,j);

        jmeastemp2 = [];
        Umeastemp2 = [];

        for ii = 1:length(U1)
            jmeastemp = useParameters(2)*(exp(useParameters(1)/(f*T)*U1(ii))-exp((useParameters(2)-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
            if jmeastemp >= useParameters(4)
                break;
            elseif jmeastemp > 0.001
                U2 = Uohm.calculate('current', jmeastemp, 'r', useParameters(3)); % Ohmic overpotential
                U3 = Ucon.calculate('current', jmeastemp, 'T', T, 'j_lim', useParameters(4)); % Concentration overpotential
                Umeastemp = Uocv.calculate()+U1(ii)+U2+U3+useParameters(5); % "Measured" voltage
                jmeastemp2 = [jmeastemp2;jmeastemp];
                Umeastemp2 = [Umeastemp2;Umeastemp];
            end
        end

        if length(jmeas) > length(jmeastemp2)
            jmeastemp2 = [jmeastemp2;nan(abs(length(jmeas)-length(jmeastemp2)),1)];
            Umeastemp2 = [Umeastemp2;nan(abs(length(Umeas)-length(Umeastemp2)),1)];
        elseif length(jmeas) < length(jmeastemp2)
            jmeas = [jmeas;nan(abs(length(jmeas)-length(jmeastemp2)),3)];
            Umeas = [Umeas;nan(abs(length(Umeas)-length(Umeastemp2)),3)];
        end

        jmeas(:,j) = jmeastemp2;
        Umeas(:,j) = Umeastemp2;
    end

    figure
    plot(jmeas,Umeas)
    title(['Test "measurement" UI, ' parameterNames(i)])
    xlabel('j')
    ylabel('U')
    legend(num2str(parameters(i,1)),num2str(parameters(i,2)),num2str(parameters(i,3)))

end
