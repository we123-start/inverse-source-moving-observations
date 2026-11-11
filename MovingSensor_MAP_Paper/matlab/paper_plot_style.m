function style = paper_plot_style(layoutTag)
% PAPER_PLOT_STYLE  Drawing parameters from the author's huatu1(3).m.
% This file ONLY defines appearance. It does not alter PDE/MAP computations.
% The original script uses figure([100 100 700 600]), automatic axes
% placement/ticks, 12-pt axes ticks, 20-pt X/Y labels and 20-pt title.
if nargin < 1, layoutTag = 'MAP_LAYOUT_single'; end %#ok<NASGU>
style.background         = [1 1 1];
style.domainColor        = [0.4 0.4 0.4];
style.domainWidth        = 1.2;
style.trueColor          = [1 0 0];
style.trueWidth          = 2;
style.estimateColor      = [0 0 1];
style.estimateWidth      = 2;
style.sensorColor        = [1 0 0];
style.sensorEdgeColor    = [0.8 0 0];
style.sensorMarkerSize   = 10;
style.sensorLineWidth    = 1.5;
style.arrowColor         = [0 0.4 1];
style.arrowWidth         = 2.5;
style.arrowHeadLength    = 0.12;
style.arrowHeadHalfWidth = 0.05;
style.arrowSamples       = 120;
style.axisRange          = [-1.1 1.1];
style.axesLineWidth      = 1;
style.gridColor          = [0.9 0.9 0.9];
style.gridAlpha          = 0.2;
style.figurePosition     = [100 100 700 600];
style.fontName           = get(0,'DefaultAxesFontName');
style.tickFont           = 12;
style.labelFont          = 20;
style.titleFont          = 20;
style.panelFont          = 12;
end
