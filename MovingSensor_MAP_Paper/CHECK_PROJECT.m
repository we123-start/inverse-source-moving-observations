function CHECK_PROJECT()
%CHECK_PROJECT Quick structural check without solving a PDE or saving files.
root = fileparts(mfilename('fullpath'));
addpath(root,'-begin');
addpath(fullfile(root,'matlab'),'-begin');
addpath(fullfile(root,'matlab','original_fem'),'-begin');
assert(exist('fmincon','file') == 2, ...
    'MATLAB Optimization Toolbox with fmincon is needed');
shapes = {'CIRCLE','KITE','FOURLEAF','PEANUT'};
for s = 1:numel(shapes)
    for e = 1:3
        fn = sprintf('%s_EXP%d.m',shapes{s},e);
        assert(exist(fullfile(root,fn),'file') == 2, ...
            'Missing entry point %s',fn);
    end
    kind=lower(shapes{s});
    f=fullfile(root,'fixtures',['V16_REPRO_NOISE_',kind,'.mat']);
    assert(exist(f,'file')==2,'Missing noise fixture: %s',f);
    Q=load(f,'mainZ','navZ');
    assert(numel(Q.mainZ)==320 && numel(Q.navZ)==6, ...
        'Invalid noise fixture: %s',f);
end
cfg=V16_CONFIG();
for s=1:numel(shapes)
    kind=lower(shapes{s});
    [c,ceq]=map_constraints(cfg.mapInitial.(kind),kind);
    assert(isempty(ceq) && all(c < -1e-8), ...
        'Initial MAP constraints fail for %s',kind);
end
% Verify the recursion itself without solving the heat equation.
for s=1:numel(shapes)
    kind=lower(shapes{s});
    m=cfg.m.(kind);speed=cfg.speed.(kind);
    startIdx=1;endIdx=cfg.mainPerWindow;
    for k=1:3
        if k<3, step=m*cfg.c1*pi;else,step=floor(m/2)*cfg.c1*pi;end
        step=step+cfg.rho0/2^(k-1);
        [nextStart,nextEnd,si,ei]=V16_NEXT_WINDOW( ...
            endIdx*cfg.dt,step/speed,cfg.dt,cfg.mainPerWindow);
        startIdx=si;endIdx=ei;
    end
    expected=266;
    if strcmp(kind,'peanut'),expected=265;end
    assert(startIdx==expected && endIdx==expected+79, ...
        'Unexpected fourth sensing window for %s',kind);
end
fprintf('PASS: 12 paper entries, 4 noise fixtures, 4 MAP initial points.\n');
fprintf('PASS: recursive fourth window = 266:345 (first 3), 265:344 (peanut).\n');
fprintf('Figures are interactive; no output folder is required.\n');
end
