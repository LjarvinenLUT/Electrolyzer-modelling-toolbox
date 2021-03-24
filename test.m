% 
close all;
clear;


%% Reversible voltage

[F,R,n_e] = get_constants;

% Measured, from CRC Handbook of Chemistry and Physics
Tmeas = [298.15;300;373.21];
deltaG = [237.141;236.839;225.160]*1e3;
Umeas = deltaG/(n_e*F);

% Modelled
T = sortrows([linspace(273,373,101) 298.15 373.21]');
U = [];
difU = [];
figure
hold on;
for i = 1:6
    U_0 = reversible(T,i);
    U = [U U_0];
    plot(T,U_0)
end
scatter(Tmeas,Umeas,'x')
hold off;
legend('1','2','3','4','5','6','tabled')
title('Reversible voltage')
xlabel('Temperature [K]')
ylabel('Voltage [V]')


% Difference of measured to modelled
ComInd = [find(T==Tmeas(1));find(T==Tmeas(2));find(T==Tmeas(3))];
UforComp = U(ComInd,:);
difU = (UforComp-Umeas);


figure
hold on
for i = 1:length(difU(1,:))
    scatter(Tmeas,difU(:,i))
end
yline(0,'--')
hold off
legend('1','2','3','4','5','6')
title('Difference to measured value')
xlabel('Temperature [K]')
ylabel('Difference in voltage [V]')

%% Water activity and vapor pressure in KOH and NaOH solutions
% From Balej's 1985 article

perT = linspace(1.5,3.5,50);
T = 1./perT.*1e3;
T = T(T<474);
perT = 1./T*1e-3;
mKOH = 2:2:18;
mNaOH = 2:2:24;

[psvKOH,aH2OKOH] = electrolyte_parameters(T,mKOH','KOH','model',1);

[psvNaOH,aH2ONaOH] = electrolyte_parameters(T,mNaOH','NaOH','model',1);

psv = water_vapor_pressure(T,'model',2);


figure
hold on;
for i = 1:numel(mKOH)
    plot(log10(psv),log10(psvKOH))
end
hold off;

figure
hold on;
for j = 1:numel(mNaOH)
    plot(log10(psv),log10(psvKOH))
end
hold off;

figure
hold on;
for i = 1:numel(mKOH)
    plot(perT,log10(aH2OKOH))
end
hold off;

figure
hold on;
for j = 1:numel(mNaOH)
    plot(perT,log10(aH2ONaOH))
end
hold off;

%% Pure water vapor pressure comparison

T = linspace(273.15,373.15,50)';

psv = NaN(length(T),3);
figure
hold on;
for i = 1:2
    psv(:,i) = water_vapor_pressure(T,'mode',i);
    plot(T,psv(:,i))
end
hold off;

