% Electrolyzer_model object testing
clear; clc; close all;
addpath("Utils");

%% Global parameters
[F,R,n_e] = get_constants();

f = R/(n_e*F);

T = 273.15+50; % Temperature
type = "PEM"; % Cell type
pH2 = 30; % bar
pO2 = 2; % bar

%% Generate test data
alpha = 0.3; %
j0 = 1e-6; % A/cm^2, exchange current density
r = 0.1; % Ohm, total resistance
j_lim = 1.5; % A/cm^2, limiting current density
Uerr = 0; % V, constant voltage error

Uocv = nernst(T,pH2,pO2,'type',type);
Uact = activation(T,'model',2);
Uohm = ohmic();
Ucon = concentration(T);

Umeas = []; % "Measured" voltage
jmeas = []; % "Measured" current based on Buttler-Volmer equation

U1 = ((0:0.001:100)*(f*T))'; % Uact/(f*T) = 0...10, Activation overpotential sweep limits

for ii = 1:length(U1)
    jmeastemp = j0*(exp(alpha/(f*T)*U1(ii))-exp((alpha-1)/(f*T)*U1(ii))); % "Measured" current based on Buttler-Volmer equation
    if jmeastemp >= j_lim
        break;
    elseif jmeastemp > 1e-5
        U2 = Uohm(r,jmeastemp); % Ohmic overpotential
        U3 = Ucon(j_lim,jmeastemp); % Concentration overpotential
        Umeastemp = Uocv+U1(ii)+U2+U3+Uerr; % "Measured" voltage
        jmeas = [jmeas;jmeastemp];
        Umeas = [Umeas;Umeastemp];
    end
end

% Take samples from the dense data vectors
N = 50; % Number of evenly taken current samples
jsamples = linspace(min(jmeas),max(jmeas),N)';
%         jsamples = linspace(min(jmeas),jL-0.01,N)'; % Excluding mass transport limitations
jmeassamp = nan(N,1);
Umeassamp = nan(N,1);
for ii = 1:N
    jdif = abs(jmeas-jsamples(ii));
    [~,ind] = min(jdif);
    jmeassamp(ii) = jmeas(ind); % Final sampled current vector
    Umeassamp(ii) = Umeas(ind); % Final sampled voltage vector
end


% Adding p*100% error to measurements
p = 0.01;
jmeassamper = jmeassamp.*(1+p*(rand(size(jmeassamp))-0.5));
Umeassamper = Umeassamp.*(1+p*(rand(size(jmeassamp))-0.5));

%% Use model class to do fitting with selected models
T = 285;
obj = electrolyzer_model("type", "alkaline", "electrolyte", "KOH") ...
    .add_overpotential(nernst(T, 1, 10, "type", "alkaline", "electrolyte", "KOH")) ...
    .add_overpotential(ohmic('type', 'alkaline')) ...
    .add_overpotential(activation(T)) ...
    .add_overpotential(concentration(T)) ...
    .combine_overpotentials() ...
    .fit_UI(Umeassamper, jmeassamper) ...
    .show_UI();

obj.fit_parameters


%%
%uniqueArguments = obj.getOverpotentialArguments()
obj.overpotential_function;
obj.overpotential_function(0.5, 1e-6, 2, 10, 4);