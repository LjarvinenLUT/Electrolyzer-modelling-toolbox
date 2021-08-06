% EMODELSTARTUP Script that adds all the necessary paths to the search
% path when called from outside the directory of the tool.

[toolpath,~,~] = fileparts(which('electrolyzerModel'));

addpath([toolpath '/Utils'])
addpath([toolpath '/Utils/mcmcstat'])
addpath([toolpath '/Utils/legendflex'])
addpath([toolpath '/Tests'])
addpath([toolpath '/Modelling'])
addpath([toolpath '/Test data'])
