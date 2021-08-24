% STARTUP Startup script that adds all the necessary paths to the search
% path.
%
% This file should be copied to any project using the electrolyzer
% modelling toolbox. This script is run automatically when matlab is opened
% in a directory containing this file. (Other commands that should be run
% when opening the project can be added in the end of this file.)

toolpath = pwd; % <--- Change the variable toolpath to a character vector 
                  %         or a string pointing to the directory where the
                  %         electrolyzer model toolbox has been saved.
addpath(toolpath) 
try
	startupEModel;
catch ME
    warningMsg = ['Path to the electrolyzer modelling toolbox not defined.'...
        'Please replace the toolpath variable with the path to the correct'...
        'directory of the toolbox.'];
    warning(warningMsg)
end