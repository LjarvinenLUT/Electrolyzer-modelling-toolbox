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
alfa = 0.5; % Electron transfer coefficient
TC = 0; % Temperature in °C
T = TC + 273.15;

Uvector = -0.1:.001:0.1;
old_a = -1111;
old_t = -1111;
legendText = [];
for i = 1:length(T)
   t = T(i);
   for j = 1:length(alfa)
       a = alfa(j);
       for k = 1:length(j0)
           j_0 = j0(k);

           % Butler-Volmer
           f = R/(n_e*F);
           BV_j = @(Uact) j_0*(exp((a*Uact)/(f*t)) - exp((a-1)*Uact/(f*t)));

           % Plot U/I curve
           if old_a ~= a || old_t ~= t
               figure("Name", "Buttler-Volmer")
               legendText = [string.empty];
           end
           hold on;
           plot(Uvector, BV_j(Uvector)); grid on;
           xlabel("Voltage (V)")
           ylabel("Current (A)")
           legendText(end+1) = "j0=" + j_0;
           legend(legendText)
           if old_a ~= a
               hold off;
           end
           old_a = a;
           old_t = t;
       end
   end
end

old_a = -1111;
old_t = -1111;
legendText = [];
for i = 1:length(T)
    t = T(i);
    for j = 1:length(alfa)
        a = alfa(j);
        for k = 1:length(j0)
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
            if  old_a ~= a || old_t ~= t
                figure("Name", "Tafel")
                legendText = [string.empty];
            end
            hold on; grid on
            plot(Uvector, real(log(Tafelvoltage)))
            plot(Uvector, real(log(Tafelvoltage1)))
            plot(Uvector, real(log(Tafelvoltage2)))
            xlabel("Voltage (V)")
            ylabel("Current (A)")
            legendText(end+1) = "j0=" + j_0;
            legend(legendText)
            title("T=" + t + " a=" + a)
            if old_a ~= a
                hold off;
            end
            old_a = a;
            old_t = t;
        end
    end
end