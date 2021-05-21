close all; clearvars; clc;
% A new test script for electrolyzerModel object

%% Constructing
                                                                                                                           
type = 'PEM';

e = electrolyzerModel('type',type)

%% Adding variables

T = 273.15 + 50; % K
pCat = 30; % bar
pAn = 2; % bar

e.setVariables('T',T,'pCat',pCat,'pAn',pAn);

e.Variables

%% Removing variables

e.setVariables('Q','Oops','p','wrong variable');
e.Variables

e.removeVariables('Q','p')
e.Variables


%% Defining potential function with defaults

e.addAllPotentials('nernst','activation','ohmic','concentration')

%% Copying the electrolyzer model

e2 = e.copy