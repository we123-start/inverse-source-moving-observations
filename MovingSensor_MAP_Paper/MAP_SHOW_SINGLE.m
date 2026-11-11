function fh = MAP_SHOW_SINGLE(kind,xi,angles,mode,windowName,idx)
%MAP_SHOW_SINGLE A single interactive MATLAB figure based on huatu1(3).m.
% This function NEVER exports, saves, closes, or crops a figure.
% Save manually using the MATLAB Figure window (File > Save As).
if nargin < 5 || isempty(windowName)
    windowName = sprintf('%s MAP',upper(kind));
end
if nargin < 6 || isempty(idx),idx = 1;end
set(0,'DefaultFigureColor','white');
set(0,'DefaultAxesFontSize',12);
set(0,'DefaultTextFontSize',12);
shift = mod(idx-1,8)*25;
% Original figure size and default axes layout from huatu1(3).m.
fh = figure('Color','white','Visible','on', ...
    'Units','pixels','Position',[100+shift,100+shift,700,600], ...
    'Name',windowName,'NumberTitle','off', ...
    'Tag','MAP_LAYOUT_single');
% Do not force axes Position or crop with a tight bounding box: huatu1(3)
% used MATLAB's default axes layout on a 700-by-600 figure.
ax = axes('Parent',fh);
plot_map_shape(ax,kind,xi(:),angles(:),mode,'');
drawnow;
end
