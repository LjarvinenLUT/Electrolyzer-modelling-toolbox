% EMODELSTARTUP Script that adds all the necessary paths to the search
% path when called from outside the directory of the tool.

[toolPath,~,~] = fileparts(which('electrolyzerModel'));

addpath([toolPath '/Utils'])
addpath([toolPath '/Utils/mcmcstat'])
addpath([toolPath '/Utils/legendflex'])
addpath([toolPath '/Tests'])
addpath([toolPath '/Modelling'])
addpath([toolPath '/Test data'])
