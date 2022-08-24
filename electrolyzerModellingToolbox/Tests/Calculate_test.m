clear; close all; clc

fprintf("\n______________OCV 1______________")
Uocv1 = nernst('PEM');
Uocv_val1 = Uocv1.calculate();

fprintf("\n______________OCV 2______________")
Uocv2 = nernst('PEM');
Uocv_val2 = Uocv2.calculate('T',300,'pCat',1);

fprintf("\n______________OCV 3______________")
Uocv3 = nernst('PEM');
Uocv_val3 = Uocv3.calculate('T',300,'pCat',1,'pAn',1);

%% Ohmic
fprintf("\n______________Ohm 1______________")
Uohm1 = ohmic();
Uohm_val1 = Uohm1.calculate();

fprintf("\n______________Ohm 2______________")
Uohm2 = ohmic();
Uohm_val2 = Uohm2.calculate('r',1);

fprintf("\n______________Ohm 3______________")
Uohm3 = ohmic();
Uohm_val3 = Uohm3.calculate('r',1,'current',1);
