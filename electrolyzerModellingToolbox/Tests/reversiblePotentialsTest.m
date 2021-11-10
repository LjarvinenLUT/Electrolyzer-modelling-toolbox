% A script to test different approximations for the reversible potential
clear; close all; clc;

T = (0:1:100) + 273.15;

f1 = figure('Name','Reversible potential approximations, values');
ax1 = gca;
hold on

f2 = figure('Name','Reversible potential approximations, difference');
ax2 = gca;
hold on


Ufun1 = nernstReversible(2);
U1 = calculate(Ufun1,'T',T);
plot(ax1,T,U1)
plot(ax2,T,U1./U1)

for i = [1 3 4 5 6]
    Ufun = nernstReversible(i);
    U = calculate(Ufun,'T',T);
    plot(ax1,T,U)
    plot(ax2,T,U./U1)
end

legend(ax1,'2','1','3','4','5','6')
legend(ax2,'2','1','3','4','5','6')