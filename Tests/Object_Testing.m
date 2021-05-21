% Electrolyzer_model object testing
clear; clc; close all;

obj = electrolyzer_model("type", "alkaline", "electrolyte", "KOH") ...
    .add_overpotential(ohmic()) ...
    .add_overpotential(ohmic()) ...
    .add_overpotential(ohmic());

uniqueArguments = obj.getOverpotentialArguments()
obj.set_overpotentials()
obj.overpotential_function(1, 1)

ohm = ohmic();
real = ohm(1,1) * 3


%% Random testing

fun = @(a) a + 1
cells = {"x"}
funString = "@(" + cells{:} + ") fun(" + cells{:} + ")"
ff = str2func(funString)
f = @(x) ff(x)
f(1) 