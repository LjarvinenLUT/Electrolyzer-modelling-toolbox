%% Function for setting figure parameters

function [fileOutputPath,FigSizeWidth,FigSizeHeight] = setFigureParameters()

if exist('FigureParameters.mat','file') == 2
    load('FigureParameters.mat');
else
    currentFolder = pwd;
    fileFormat = '.pdf';
    fileOutputPath = [pwd '\Figures\'];
    
    % Figure parameters
    FontSize = 8;
    FontName = 'Times New Roman';
    LineWidth = 1;
    Units = 'centimeters';
    Location = 'SouthWest';
    
    % Monitor size
    FigSizeWidth = 9;%8.5 %IEEE 3.5in ~ 8.89cm Elsevier 9cm
    %https://journals.ieeeauthorcenter.ieee.org/create-your-ieee-journal-article/create-graphics-for-your-article/accepted-graphics-file-formats/
    FigSizeHeight = 7;     % One plot
    FigPosLeft = 0;
    FigPosBottom = 19;
    GridLineStyle = ':';
    
    % LUT Colors
    % Red R237 G23 B77
    LUTred = [237/255 23/255 77/255];
    % Orange R255 G102 B0
    LUTorange = [255/255 102/255 0];
    % Green R64 G213 B32
    LUTgreen = [64/255 213/255 32/255];
    % Black R0 G0 B0
    LUTblack = [0 0 0];
    % Blue R0 G102 B255
    LUTblue = [0 102/255 255/255];
    
    LUTmap = [LUTorange;LUTgreen;LUTred;LUTblue;LUTblack];
    
    % addpath(genpath('\ghostscript files'));
    
    
    marks = {'o','+','*','x','s','d','^','v','p','h'};
    fc = {'flat','none','none','none','flat','flat','flat','flat','flat','flat'}; % Face color
    ec = {'none','flat','flat','flat','none','none','none','none','none','none'}; % Edge color
    
    
    save('FigureParameters.mat',...
        'fileFormat',...
        'fileOutputPath',...
        'FontName',...
        'FontSize',...
        'FigSizeWidth',...
        'FigSizeHeight',...
        'FigPosLeft',...
        'FigPosBottom',...
        'GridLineStyle',...
        'Units',...
        'LUTred',...
        'LUTorange',...
        'LUTgreen',...
        'LUTblack',...
        'LUTblue',...
        'LUTmap',...
        'marks',...
        'fc',...
        'ec');
end

set(groot,'defaultAxesLineStyleOrder','-');
% set(groot,'defaultAxesColorOrder',LUTmap);
set(groot,'defaultFigureUnits',Units); 
set(groot,'defaultFigurePosition',[FigPosLeft FigPosBottom FigSizeWidth FigSizeHeight]); 
set(groot,'defaultAxesFontName', FontName);
set(groot,'defaultAxesFontSize', FontSize);
set(groot,'defaultAxesGridLineStyle', GridLineStyle);
set(groot,'defaultFigureColor', 'w');
set(groot,'defaultAxesFontName', FontName);
set(groot,'defaultAxesFontSize', FontSize);
%set(groot,'LegendInterpreter', 'latex');
set(groot,'defaultTextInterpreter','latex');
set(groot,'defaultLegendInterpreter','latex');
set(groot,'defaultAxesTickLabelInterpreter','latex');

end