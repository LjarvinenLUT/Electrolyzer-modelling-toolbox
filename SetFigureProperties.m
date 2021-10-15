function SetFigureProperties()
%SetFigureProperties Define the default figure properties suitable for
%publications

%% Set units
Units = 'centimeters';
set(groot,'defaultFigureUnits',Units);

%% Monitor size
set(0,'units',Units)
screenSize = get(0,'ScreenSize'); % Screen size [~ ~ width height]

% Figure size and position
FigSizeWidth = 9;%8.5 %IEEE 3.5in ~ 8.89cm Elsevier 9cm
% %https://journals.ieeeauthorcenter.ieee.org/create-your-ieee-journal-article/create-graphics-for-your-article/accepted-graphics-file-formats/
FigSizeHeight = 7;     % One plot
FigPosLeft = 0;
FigPosBottom = 0;

set(groot,'defaultFigurePosition',[FigPosLeft FigPosBottom FigSizeWidth FigSizeHeight]); 

f = figure;
outPos = get(f,'OuterPosition');
close(f);
FigPosLeft = FigPosLeft-outPos(1);
FigPosBottom = screenSize(4)-outPos(4);
set(groot,'defaultFigurePosition',[FigPosLeft FigPosBottom FigSizeWidth FigSizeHeight]); 

%% Figure parameters
FontSize = 8;
FontName = 'Times New Roman';
LineWidth = 1;
Location = 'SouthWest';
GridLineStyle = ':';

set(groot,'defaultAxesFontName', FontName);
set(groot,'defaultAxesFontSize', FontSize);
set(groot,'defaultAxesGridLineStyle', GridLineStyle);
set(groot,'defaultFigureColor', 'w');
set(groot,'defaultAxesFontName', FontName);
set(groot,'defaultAxesFontSize', FontSize);
%set(groot,'LegendInterpreter', 'latex');
%set(groot,'defaultTextInterpreter','latex');
%set(groot,'defaultLegendInterpreter','latex');
%set(groot,'defaultAxesTickLabelInterpreter','latex');
end

