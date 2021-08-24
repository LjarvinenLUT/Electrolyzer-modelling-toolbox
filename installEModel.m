% INSTALLEMODEL Script that adds all the necessary paths to the search
%   path when called from outside the directory of the tool and saves them
%   to the default search path.

[toolPath,~,~] = fileparts(which('electrolyzerModel'));

addpath([toolPath '/Utils'])
addpath([toolPath '/Utils/mcmcstat'])
addpath([toolPath '/Utils/legendflex'])
addpath([toolPath '/Tests'])
addpath([toolPath '/Modelling'])
addpath([toolPath '/TestData'])

savepath;