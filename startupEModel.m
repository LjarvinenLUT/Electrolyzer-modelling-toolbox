% STARTUPEMODEL Script that adds all the necessary paths to the search
%   path when called from outside the directory of the tool.

[toolpath,~,~] = fileparts(which('electrolyzerModel'));

addpath(genpath([toolpath '/Utils']))
addpath([toolpath '/Tests'])
addpath([toolpath '/Modelling'])
addpath([toolpath '/TestData'])
addpath(genpath([toolpath '/Examples']))
addpath([toolpath '/doc'])
