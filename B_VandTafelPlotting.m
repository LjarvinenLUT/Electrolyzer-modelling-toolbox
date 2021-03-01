close all;
clear;

% Globals dont seem to work for some reason so everything is defined here
% TODO: FIX
N_A = 6.02214076e23; % 1/mol, Avogadro constant
k_B = 1.380649e-23; % J/K, Boltzmann constant
Q_e = 1.602176634e-19; % C, Elemental charge
F = Q_e*N_A; % C/mol, Faraday constant
R = k_B*N_A; % J/K/mol, Universal gas constant
n_e = 2; % Number of electrons transferred in one reaction

j0 = 1e-3:.001:5e-3; % Exchange current density
alfa = 0.5:0.05:0.7; % Electron transfer coefficient
TC = 0:20:80; % Temperature in °C
T = TC + 273.15;

Uvector = -0.1:.001:0.1;
old_a = -1111;
old_t = -1111;

% BV Plotting changing temperature (T)
figure("Name", "Buttler-Volmer T")
legendText = [string.empty];
hold on;
for i = 1:length(T)
   t = T(i);
   a = alfa(1);
   j_0 = j0(1);
   
   % Butler-Volmer
   f = R/(n_e*F);
   BV_j = @(Uact) j_0*(exp((a*Uact)/(f*t)) - exp((a-1)*Uact/(f*t)));
   
   % Plot U/I curve
   plot(Uvector, BV_j(Uvector)); grid on;
   xlabel("Voltage (V)")
   ylabel("Current (A)")
   title("a=" + a + " j0=" + j_0)
   legendText(end+1) = "T=" + t;
   legend(legendText)
end
hold off;

% BV Plotting changing charge transfer coefficient (alfa)
figure("Name", "Buttler-Volmer alfa")
legendText = [string.empty];
hold on;
for j = 1:length(alfa)
    t = T(1);
    a = alfa(j);
    j_0 = j0(1);
    
    % Butler-Volmer
    f = R/(n_e*F);
    BV_j = @(Uact) j_0*(exp((a*Uact)/(f*t)) - exp((a-1)*Uact/(f*t)));
    
    % Plot U/I curve
    plot(Uvector, BV_j(Uvector)); grid on;
    xlabel("Voltage (V)")
    ylabel("Current (A)")
    title("T=" + t + " j0=" + j_0)
    legendText(end+1) = "a=" + a;
    legend(legendText)
end
hold off

% BV Plotting changing exchange current density (j0)
figure("Name", "Buttler-Volmer j0")
legendText = [string.empty];
hold on;
for k = 1:length(j0)
    t = T(1);
    a = alfa(1);
    j_0 = j0(k);
    
    % Butler-Volmer
    f = R/(n_e*F);
    BV_j = @(Uact) j_0*(exp((a*Uact)/(f*t)) - exp((a-1)*Uact/(f*t)));
    
    % Plot U/I curve
    plot(Uvector, BV_j(Uvector)); grid on;
    xlabel("Voltage (V)")
    ylabel("Current (A)")
    title("T=" + t + " a=" + a)
    legendText(end+1) = "j0=" + j_0;
    legend(legendText)
end
hold off;

% Tafel Plotting changing exchange current density (j0)
figure("Name", "Tafel j0")
legendText = [string.empty];
hold on; grid on
for k = 1:length(j0)
    t = T(1);
    a = alfa(1);
    j_0 = j0(k);
    % Tafel equation for both sides
    Tafel_jp = @(Uact) j_0*exp((n_e*(1 - a)*F.*Uact)/(R*t));
    Tafel_jn = @(Uact) -j_0*exp((n_e*-a.*F*Uact)/(R*t));
    
    Tafelvoltage1 = [];
    Tafelvoltage2 = [];
    for i = 1:length(Uvector)
        U = Uvector(i);
        Tafelvoltage1(end+1) = Tafel_jp(U);
        Tafelvoltage2(end+1) = Tafel_jn(U);
    end
    
    % Calculate combined voltage
    Tafelvoltage = Tafelvoltage1 + Tafelvoltage2;
    
    % Plot U/I curve
    plot(Uvector, real(log(Tafelvoltage)))
    plot(Uvector, real(log(Tafelvoltage1)))
    plot(Uvector, real(log(Tafelvoltage2)))
    xlabel("Voltage (V)")
    ylabel("Current (A)")
    legendText(end+1) = "j0=" + j_0;
    legend(legendText)
    title("T=" + t + " a=" + a)
end
hold off;

% Tafel Plotting changing temperature (T)
figure("Name", "Tafel T")
legendText = [string.empty];
hold on; grid on
for i = 1:length(T)
    t = T(i);
    a = alfa(1);
    j_0 = j0(1);
    % Tafel equation for both sides
    Tafel_jp = @(Uact) j_0*exp((n_e*(1 - a)*F.*Uact)/(R*t));
    Tafel_jn = @(Uact) -j_0*exp((n_e*-a.*F*Uact)/(R*t));
    
    Tafelvoltage1 = [];
    Tafelvoltage2 = [];
    for i = 1:length(Uvector)
        U = Uvector(i);
        Tafelvoltage1(end+1) = Tafel_jp(U);
        Tafelvoltage2(end+1) = Tafel_jn(U);
    end
    
    % Calculate combined voltage
    Tafelvoltage = Tafelvoltage1 + Tafelvoltage2;
    
    % Plot U/I curve
    plot(Uvector, real(log(Tafelvoltage)))
    plot(Uvector, real(log(Tafelvoltage1)))
    plot(Uvector, real(log(Tafelvoltage2)))
    xlabel("Voltage (V)")
    ylabel("Current (A)")
    legendText(end+1) = "T=" + t;
    legend(legendText)
    title("a=" + a + " j0=" + j_0)
end
hold off;

% Tafel Plotting changing charge transfer coefficient (alfa)
figure("Name", "Tafel a")
legendText = [string.empty];
hold on; grid on
for j = 1:length(alfa)
    t = T(1);
    a = alfa(j);
    j_0 = j0(1);
    % Tafel equation for both sides
    Tafel_jp = @(Uact) j_0*exp((n_e*(1 - a)*F.*Uact)/(R*t));
    Tafel_jn = @(Uact) -j_0*exp((n_e*-a.*F*Uact)/(R*t));
    
    Tafelvoltage1 = [];
    Tafelvoltage2 = [];
    for i = 1:length(Uvector)
        U = Uvector(i);
        Tafelvoltage1(end+1) = Tafel_jp(U);
        Tafelvoltage2(end+1) = Tafel_jn(U);
    end
    
    % Calculate combined voltage
    Tafelvoltage = Tafelvoltage1 + Tafelvoltage2;
    
    % Plot U/I curve
    plot(Uvector, real(log(Tafelvoltage)))
    plot(Uvector, real(log(Tafelvoltage1)))
    plot(Uvector, real(log(Tafelvoltage2)))
    xlabel("Voltage (V)")
    ylabel("Current (A)")
    legendText(end+1) = "a=" + a;
    legend(legendText)
    title("T=" + t + " j0=" + j_0)
end
hold off;