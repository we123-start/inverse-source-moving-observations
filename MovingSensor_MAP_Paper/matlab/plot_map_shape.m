function info = plot_map_shape(ax,kind,xi,angles,mode,panelLabel)
% Geometry plot styled to match huatu(2), huatu1(2), huatu2(2).
% All MAP xi values and sensor angles are supplied by the experiment data.
%
% Modes:
%   stage: stage 1 has one sensor; stages 2--4 show ONLY the latest two.
%   policy: single fixed sensor, no movement arrow.
%   comparison_moving: two stationary sensors and one directed blue arc.
%   adaptive_final_compare: show all FOUR adaptive sites, with NO arrows.
%       This is the final, publication-sized Exp-1 comparison figure.
if nargin<6, panelLabel=''; end
assert(ishandle(ax),'Invalid axes handle');
assert(any(strcmp(mode,{'stage','policy','comparison_moving', ...
    'adaptive_final_compare'})), ...
    'Unknown plot mode');
angles=angles(:);
assert(~isempty(angles) && all(isfinite(angles)),'Invalid sensor angle input');

shapes={'circle','kite','fourleaf','peanut'};
titles={'Circle','Kite','4-Leaf','Peanut'};
idx=find(strcmpi(kind,shapes),1);
assert(~isempty(idx),'Unknown geometry %s',kind);
fig=ancestor(ax,'figure');
style=paper_plot_style(get(fig,'Tag'));

% Correct geometric parametric boundaries for each individual problem.
% No hand-crafted parameter vectors are inserted here.
theta=linspace(0,2*pi,400);
truth=reviewer2_shapes('truth',kind);
actual=reviewer2_shapes('boundary',kind,truth,theta);
estimate=reviewer2_shapes('boundary',kind,xi(:),theta);

cla(ax);
hold(ax,'on');
set(ax,'Color',style.background,'FontName',style.fontName, ...
    'FontSize',style.tickFont,'LineWidth',style.axesLineWidth, ...
    'Box','on','XGrid','off','YGrid','off');

% Original huatu colors / widths, kept exactly.
plot(ax,cos(theta),sin(theta),'-', ...
    'Color',style.domainColor,'LineWidth',style.domainWidth, ...
    'Tag','MAP_DomainBoundary');
plot(ax,actual(1,:),actual(2,:),'r--', ...
    'LineWidth',style.trueWidth,'Tag','MAP_TrueBoundary');
plot(ax,estimate(1,:),estimate(2,:),'b:', ...
    'LineWidth',style.estimateWidth,'Tag','MAP_EstimatedBoundary');

% There are 80 data at each stationary site, but the location trajectory
% itself comes exclusively from the stored original time/angle arrays.
assert(mod(numel(angles),80)==0, ...
    'Expected exactly 80 main observations per stationary window');
stationary=angles(1:80:end);
arrowPair=[];
if strcmp(mode,'stage')
    stage=numel(stationary);
    assert(stage>=1 && stage<=4,'Stage must be 1 to 4');
    if stage==1
        selected=stationary(1);
    else
        selected=stationary([stage-1 stage]);
        arrowPair=[1 2];
    end
elseif strcmp(mode,'adaptive_final_compare')
    % The comparison uses the Stage-4 estimate and ALL 320 main readings.
    % Display all four locations; do not draw any movement arrows.
    assert(numel(angles)==320 && numel(stationary)==4, ...
        'Final adaptive figure requires four 80-reading windows');
    selected=stationary(:).';
    for jj=1:4
        other=selected([1:jj-1 jj+1:4]);
        assert(all(abs(atan2(sin(selected(jj)-other), ...
            cos(selected(jj)-other)))>1e-7), ...
            'Final adaptive observation sites must be distinct');
    end
    % Keep arrowPair empty: four markers and no circular movement arc.
else
    selected=local_distinct_angles(stationary,1e-7);
    if strcmp(mode,'comparison_moving')
        assert(numel(selected)==2, ...
            'Two stationary sites are required for a local-moving figure');
        arrowPair=[1 2];
    else
        assert(numel(selected)==1, ...
            'A fixed-location figure requires one stationary site');
    end
end
showArrow=~isempty(arrowPair);

% Draw arcs (when requested) BEFORE red markers, so sites stay visible.
% The final four-site comparison deliberately has no arrows. Stage and
% local two-site movement figures keep their original arcs.
if showArrow
    local_arc_arrow(ax,selected(arrowPair(1)),selected(arrowPair(2)),style);
end
for k=1:numel(selected)
    plot(ax,cos(selected(k)),sin(selected(k)),'o', ...
        'MarkerSize',style.sensorMarkerSize, ...
        'LineWidth',style.sensorLineWidth, ...
        'MarkerFaceColor',style.sensorColor, ...
        'MarkerEdgeColor',style.sensorEdgeColor, ...
        'Tag','MAP_ObservationPoint');
end

axis(ax,'equal');
xlim(ax,style.axisRange);ylim(ax,style.axisRange);
set(ax,'GridColor',style.gridColor,'GridAlpha',style.gridAlpha);
grid(ax,'off');box(ax,'on');
xlabel(ax,'X-axis','FontSize',style.labelFont,'FontName',style.fontName);
ylabel(ax,'Y-axis','FontSize',style.labelFont,'FontName',style.fontName);
title(ax,['Actual vs Estimated ',titles{idx}], ...
    'FontSize',style.titleFont,'FontName',style.fontName);
if ~isempty(panelLabel)
    text(ax,0.03,0.965,panelLabel,'Units','normalized', ...
        'HorizontalAlignment','left','VerticalAlignment','top', ...
        'FontSize',style.panelFont,'FontName',style.fontName, ...
        'Interpreter','none','Color',[0.1 0.1 0.1], ...
        'Tag','MAP_PanelLabel');
end
hold(ax,'off');

% Actual plot-object audits prevent a picture silently missing its arrow.
points=findobj(ax,'Tag','MAP_ObservationPoint');
arcs=findobj(ax,'Tag','MAP_StageMovement');
heads=findobj(ax,'Tag','MAP_Arrowhead');
assert(numel(points)==numel(selected),'Missing observation markers');
if showArrow
    assert(numel(arcs)==1 && numel(heads)==2, ...
        'Expected one blue movement arc and two open arrowhead lines');
else
    assert(isempty(arcs) && isempty(heads), ...
        'This non-moving plot mode must not contain any arrows');
end
info=struct('nMarkers',numel(points), ...
    'selectedAngles',selected(:),'hasArrow',logical(showArrow), ...
    'nArcs',numel(arcs),'nArrowheadLines',numel(heads), ...
    'arrowPair',arrowPair);
if showArrow
    info.arrowStartAngle=selected(arrowPair(1));
    info.arrowEndAngle=selected(arrowPair(2));
else
    info.arrowStartAngle=NaN;
    info.arrowEndAngle=NaN;
end
end

function uniq=local_distinct_angles(angles,tol)
uniq=[];
for k=1:numel(angles)
    if isempty(uniq) || all(abs(atan2(sin(angles(k)-uniq), ...
                                      cos(angles(k)-uniq)))>tol)
        uniq(end+1)=angles(k); %#ok<AGROW>
    end
end
end

function local_arc_arrow(ax,startAngle,endAngle,style)
% EXACT reference style: 120 arc points, bright blue dashed line width
% 2.5, classic open arrowhead length 0.12 and half-width 0.05.
% Unlike the reference examples, keep the TRUE signed movement direction.
delta=atan2(sin(endAngle-startAngle),cos(endAngle-startAngle));
assert(abs(delta)>1e-8,'The two observation positions coincide');
a=linspace(startAngle,startAngle+delta,style.arrowSamples);
arcX=cos(a);arcY=sin(a);
plot(ax,arcX,arcY,'--','Color',style.arrowColor, ...
    'LineWidth',style.arrowWidth,'Tag','MAP_StageMovement');

% Reference huatu1/2: use end-3 of arc to form the tangent-like direction.
base=[arcX(end-3);arcY(end-3)];
tip=[arcX(end);arcY(end)];
v=tip-base;
v=v/norm(v);
perp=[-v(2);v(1)];
L=style.arrowHeadLength;
W=style.arrowHeadHalfWidth;
left=tip-L*v+W*perp;
right=tip-L*v-W*perp;
plot(ax,[left(1) tip(1)],[left(2) tip(2)],'-', ...
    'Color',style.arrowColor,'LineWidth',style.arrowWidth, ...
    'Tag','MAP_Arrowhead');
plot(ax,[right(1) tip(1)],[right(2) tip(2)],'-', ...
    'Color',style.arrowColor,'LineWidth',style.arrowWidth, ...
    'Tag','MAP_Arrowhead');
end
