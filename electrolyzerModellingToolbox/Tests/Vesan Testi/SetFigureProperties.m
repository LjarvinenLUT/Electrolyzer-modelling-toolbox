function SetFigureProperties()
%SetFigureProperties Define the default figure properties suitable for
%publications

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

%set(groot,'defaultAxesColorOrder',LUTmap,...
%      'defaultAxesLineStyleOrder','-x|-o|:');
set(groot,'defaultFigureUnits',Units); 
set(groot,'defaultFigurePosition',[FigPosLeft FigPosBottom FigSizeWidth FigSizeHeight]); 
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

