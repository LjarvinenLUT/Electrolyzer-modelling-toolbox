close all; clearvars; clc;
% Test script for all the overpotential functions

%% Reversible potential

T = 273.15+(0:80);

for i = 1:6
    Urev = reversible(T,i)
end

%% Nernst equation pressure correction

p1 = 2;
p2 = 30;
m = 7.5;
%% PEM
Urev = nernstPressureCorrection(T,p1,p2,'type','pem')
%% Alkali
Urev = nernstPressureCorrection(T,p1,m,'type','alkaline')

Urev = nernstPressureCorrection(T,p1,m,'type','alkali')

%% Nernst equation

%% PEM
Uocv = nernst(T,p1,p2,'type','pem')

%% Alkali
Uocv = nernst(T,p1,m,'type','alkali')

%% Activation overpotential

Uact = activation('model',1)

%%
Uact = activation('model',2)

%%
Uact = activation('model',3)

%% Ohmic overpotential

Uohm = ohmic;

%% 
Uohm = ohmic('type','pem','resistanceModel',2,'lambda',15,'delta',2,'Temperature',T)