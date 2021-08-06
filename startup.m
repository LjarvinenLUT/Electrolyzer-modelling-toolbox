% STARTUP Startup script that adds all the necessary paths to the search
% path

toolPath = pwd;
addpath(toolPath)
try
	startupEModel;
catch ME
    warningMsg = ['Path to the electrolyzer modelling toolbox not defined.'...
        'Please replace the toolpath variable with the path to the correct'...
        'directory of the toolbox.'];
    warning(warningMsg)
end