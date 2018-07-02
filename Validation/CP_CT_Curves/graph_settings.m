function []=graph_settings
% The following changes the defualt matlab setting for the current matlab
% session. If you logout and in again they reset to defaults.
% I have set this up in the hope that it is Ignazio friendly and designed
% to make nice plots in a report/paper.
% I would reccommend downloading the export_fig function and saving your
% plots directly with the command: export_fig -transparent -nocrop path/filname.pdf for the best results.
 
% Set the colours, line porperties & font size of multi-line plots to IMV pref
% colour RGB values (7 colours)

H=[ 0         0.4470    0.7410
    0.8500    0.3250    0.0980
    0         0.5       0
    0.4940    0.1840    0.5560
    0.4660    0.6740    0.1880
    0.3010    0.7450    0.9330
    0.6350    0.0780    0.1840];
set(0,'DefaultAxesColorOrder',H)
style={'s';'o';'x';'p';'*';'v';'d';'^';'h';'+';'>';'<';'.'};
% line width is higher since most often plots are scaled down in papers
set(0,'DefaultLineLineWidth',1.5,'DefaultAxesFontSize',20,...
    'DefaultTextFontSize',16)
set(0,'defaultfigurecolor',[1 1 1]) % set frame to white
 
% Allows you to use latex formatting on figures
set(0,'defaulttextinterpreter','latex')
set(0,'DefaultLegendInterpreter','latex')
set(0,'defaultAxesTickLabelInterpreter','latex');
 
 
set(0,'DefaultAxesPosition',[0.17 0.15 0.66 0.76]) % controls size of inner plot box, see note
set(0,'DefaultFigurePosition', [100, 100, 550, 500]) % controls size of outer plot, inc whitespace, see note
%       Note: *Position is [xcord ycord xlength ylength]
%       DefaultAxesPosition xcord is 0.17*500 as the plot window is 550 wide.
%       Plot window set to be larger than normal to allow for equations on
%       axis labels. When you save figures they are alligned to the inner
%       plot box, not outer plot window.
set(0,'DefaultAxesXGrid','on','DefaultAxesYGrid','on','DefaultAxesZGrid','on') % turn grid on
set(0,'DefaultAxesBox','on') % turn box on

end