clear;
close all;
clc;

%% Global parameters
[F,R,n_e] = get_constants();

f = R/(n_e*F);


%% Creating function handle for fit

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pH2 = 30; % bar
pO2 = 2; % bar
Uocv = nernst(T,pH2,pO2,'type',type);
Uact = activation(T);
Uohm = ohmic();
Ucon = concentration(T);

Uforfit = @(j0,alpha,r,jL,Uerr,Current) Uocv + Uact(j0,alpha,Current) + Uohm(r,Current) + Ucon(jL,Current) + Uerr;

%% Creating test data

%             U1 = ((0.01:0.01:32)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential
%             Tset = ones(size(U1))*T; % T vector
%             alpha = 0.4;
%             j0 = 1.3e-6;
%             x = j0*(exp(alpha./(f*Tset).*U1)-exp((alpha-1)./(f*Tset).*U1)); % "Measured" current based on Buttler-Volmer equation
%             r = 0.3;
%             U2 = r*x; % Ohmic overpotential
%             jL = 1.5;
%             U3 = real(Ucon(jL,x));
%
%             U4 = ones(size(U1))*Uerr; % Constant potential error
%             U0 = ones(size(U1))*Uocv; % Open circuit voltage
%             z = U0+U1+U2+U3+U4; % "Measured" voltage

alpha = 0.2719; %
j0 = 1.3533e-6; % A/cm^2, exchange current density
r = 0.3631; % Ohm, total resistance
jL = 1.4092; % A/cm^2, limiting current density
Uerr = 0; % V, constant voltage error

Umeas = []; % "Measured" voltage
jmeas = []; % "Measured" current based on Buttler-Volmer equation

U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits

for ii = 1:length(U1)
    jmeastemp = j0*(exp(alpha/(f*T)*U1(ii))-exp((alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
    if jmeastemp >= jL
        break;
    elseif jmeastemp > 0.001
        U2 = Uohm(r,jmeastemp); % Ohmic overpotential
        U3 = Ucon(jL,jmeastemp); % Concentration overpotential
        Umeastemp = Uocv+U1(ii)+U2+U3+Uerr; % "Measured" voltage
        jmeas = [jmeas;jmeastemp];
        Umeas = [Umeas;Umeastemp];
    end
end

% Take samples from the dense data vectors
N = 600; % Number of evenly taken current samples
jsamples = linspace(min(jmeas),max(jmeas),N)';
jmeassamp = nan(N,1);
Umeassamp = nan(N,1);
for ii = 1:N
    jdif = abs(jmeas-jsamples(ii));
    [~,ind] = min(jdif);
    jmeassamp(ii) = jmeas(ind); % Final sampled current vector
    Umeassamp(ii) = Umeas(ind); % Final sampled voltage vector
end


% Adding p*100% error to measurements
p = 0;
jmeassamp = jmeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
Umeassamp = Umeassamp.*(1+p*(rand(size(jmeassamp))-0.5));


figure
scatter(jmeassamp,Umeassamp)
hold on
plot(jmeas,Umeas)
title('Test "measurement" UI')
xlabel('j')
ylabel('U')

%% Fit

fit_param = fit_UI(Uforfit,Umeassamp,jmeassamp);

j0fit = fit_param(1);
alphafit = fit_param(2);
rfit = fit_param(3);
jLfit = fit_param(4);
Uerrfit = fit_param(5);
% j0fit = j0;
% alphafit = alpha;
% rfit = r;
% jLfit = jL;
% Uerrfit = Uerr;

Ufit = Uforfit(j0fit,alphafit,rfit,jLfit,Uerrfit,jmeas);

MSE = mean((Ufit-Umeas).^2); % Mean squares error

% Plotting

figure
hold on;
scatter(jmeassamp,Umeassamp)
plot(jmeas,Ufit)
xlabel("j (A)")
ylabel("U (V)")
legend("Data", "Fit", "Location", "Best")
hold off;


