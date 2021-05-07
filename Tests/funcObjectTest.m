close all; clear; clc;

%% Func object test script
addpath('Utils')

Coefficients1 = struct('a',3,'b',12);
Variables1 = struct('I',1:30,'V',linspace(0,2,30));
Constants1 = struct('twoPi',2*pi);
Workspace1 = struct('Coefficients',Coefficients1,'Variables',Variables1,'Constants',Constants1);

funcHandle1 = @(Workspace) Workspace.Constants.twoPi*Workspace.Coefficients.a^2 + Workspace.Variables.I.*Workspace.Variables.V*Workspace.Coefficients.b;

func1 = func(funcHandle1,Workspace1);

result = calculate(func1)

fprintf(['\n' func1.getEquation '\n'])

%% Combine two func objects

Coefficients2 = struct('a',3,'b',12,'c',60);
Variables2 = struct('T',285,'p',20);
Workspace2 = struct('Coefficients',Coefficients2,'Variables',Variables2);

funcHandle2 = @(Workspace) Workspace.Coefficients.a^(1/2) + Workspace.Variables.T./Workspace.Variables.p*Workspace.Coefficients.c;

func2 = func(funcHandle2,Workspace2);

func3 = addFuncs(func1,func2)

fprintf(['\n' func3.getEquation '\n'])