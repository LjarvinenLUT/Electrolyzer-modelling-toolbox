close all; clearvars; clc;
% A new test script for electrolyzerModel object

%% Constructing
                                                                                                                           
type = 'PEM';

e = electrolyzerModel('type',type)

%% Adding variables

T = 273.15 + 50; % K
pCat = 30; % bar
pAn = 2; % bar

Workspace = struct('Variables',struct('T',T,'pCat',pCat,'pAn',pAn),'Constants',struct('Q','Oops','p','wrong variable'));

func.isWorkspace(Workspace)

e.setParams(Workspace);

e.viewWorkspace;

%% Removing variables

e.removeParams('Q','p')
e.viewWorkspace;

%% Replacing variables

e.replaceParams('T',T+30)
e.viewWorkspace;

%% Defining potential function with defaults

e.addPotentials('nernst','activation','ohmic','concentration')

disp(e.funcStorage)

%% Copying the electrolyzer model

e2 = e.copy
e2.clearPotentials

%% Adding potential with 'names'

e2.addPotentials('nernst','act',func.createEmpty,'names',{'ocv','Cat'})

disp(e2.funcStorage)




