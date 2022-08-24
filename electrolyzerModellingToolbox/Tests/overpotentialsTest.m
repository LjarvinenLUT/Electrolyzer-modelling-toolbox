close all; clearvars; clc;
% Test script for all the overpotential functions

%% Reversible potential

T = 273.15+(0:80);

for i = 1:6
    Urev = nernstReversible(1)
end

%% Nernst equation pressure correction

p1 = 2;
p2 = 30;
m = 7.5;
% %% PEM
Ucor = nernstPressureCorrection("PEM")
% %% Alkali
% Urev = nernstPressureCorrection(T,p1,m,'type','alkaline')
% 
% Urev = nernstPressureCorrection(T,p1,m,'type','alkali')

%% Nernst equation

%% PEM
Uocv = nernst("PEM")

% %% Alkali
% Uocv = nernst(T,p1,m,'type','alkali')

%% Activation overpotential

% Uact = activation('model',)
% 
% %%
Uact = activation('model',2)
% 
% %%
% Uact = activation('model',3)

%% Ohmic overpotential

Uohm = ohmic

% %% 
% Uohm = ohmic('type','pem','resistanceModel',2,'lambda',15,'delta',2,'Temperature',T)

%% Concentration overpotential
Ucon = concentration()


%% Combine overpotentials

U = addFuncs(addFuncs(Uocv,Uact),addFuncs(Uohm,Ucon))

%U.Workspace.Variables.current = linspace(0.1,10,length(T));

%% Destructurize function handle for fitting
[UFuncHandle,coeffNames,problemNames,problem] = U.destructurize('current')


%% Calculate
alpha = 0.4;
j0 = 1e-5;
j_lim = 1.5;
r = 5;
current = linspace(0.01,j_lim-0.01,81);

tic
result = U.calculate('alpha',alpha,'j0',j0,'j_lim',j_lim,'r',r,'current',current);
toc

