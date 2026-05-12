function MAE6760startup(font_size);
%
% define colors for plotting
global MCcolors; %define colors as global with access
MCcolors.red=[200,0,0]/255;
MCcolors.blue=[4,51,255]/255;
MCcolors.purple=[147,23,255]/255;
MCcolors.green=[0,160,0]/255;
MCcolors.orange=[253,128,8]/255;
MCcolors.mag=[255,64,255]/255;
MCcolors.cyan=[0,230,255]/255;
%
% define standard figure positioning and size
set(groot,'DefaultFigureUnits','pixels');
set(groot,'DefaultFigurePosition',[100 100 1600 600]);
set(groot,'DefaultFigureWindowStyle','normal');  % Important
set(groot,'DefaultAxesFontSize',16);
set(groot,'DefaultAxesFontWeight','bold');
set(groot,'DefaultLineLineWidth',2);
%
end