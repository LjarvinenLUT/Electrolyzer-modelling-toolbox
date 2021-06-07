close all; clear; clc;

%% Func object test script

Coefficients1 = struct('a',3,'b',12);
Variables1 = struct('I',1:30,'V',linspace(0,2,30));
Constants1 = struct('twoPi',2*pi);
Workspace1 = struct('Coefficients',Coefficients1,'Variables',Variables1,'Constants',Constants1);

funcHandle1 = @(Workspace) Workspace.Constants.twoPi*Workspace.Coefficients.a^2 + Workspace.Variables.I.*Workspace.Variables.V*Workspace.Coefficients.b;

func1 = func(funcHandle1,Workspace1);

result = calculate(func1)

fprintf(['\n' func1.equation '\n'])

%% Combine two func objects

Coefficients2 = struct('a',3,'b',12,'c',60);
Variables2 = struct('T',285,'p',20);
Workspace2 = struct('Coefficients',Coefficients2,'Variables',Variables2);

funcHandle2 = @(Workspace) Workspace.Coefficients.a^(1/2) + Workspace.Variables.T./Workspace.Variables.p*Workspace.Coefficients.c;

func2 = func.createEmpty;

func2.setFuncHandle(funcHandle2);

func2.setParams(Workspace2);

disp(func2.viewWorkspace);

func2.removeParams('b');

disp(func2.viewWorkspace);

func2.replaceParams('a',15);

disp(func2.viewWorkspace);

func3 = func.add(func1,func2)

fprintf(['\n' func3.equation '\n'])

disp(func3.viewWorkspace);