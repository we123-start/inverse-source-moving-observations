function out = FOURLEAF_EXP1(noiseCase,rep)
%FOURLEAF_EXP1 Independent Experiment 1 for the fourleaf source.
% Original Q2 polar FEM + backward Euler + constrained SQP MAP.
% Usage: FOURLEAF_EXP1                     (original fixed noisy realization)
%        FOURLEAF_EXP1('zero',1)           (zero-noise experiment)
% Each call recomputes its OWN measurements/solutions from scratch.
% All figures stay OPEN; images and MAT results are NEVER auto-saved.
if nargin < 1 || isempty(noiseCase), noiseCase = 'noisy'; end
if nargin < 2 || isempty(rep), rep = 1; end
S = MAP_PREPARE_CASE('fourleaf',noiseCase,rep);
cfg = S.cfg;
w0 = S.w0;
R = struct();
F = cell(1,4);
% Four sequential MAP reconstructions with 80,160,240,320 main data.
stages = cell(1,4);
for k = 1:4
    d = MAP_STAGE_DATA(S.data.adaptive,k,cfg);
    stages{k} = MAP_SOLVE_MEMORY(S.invModel,'fourleaf',d,w0,cfg);
    MAP_PRINT_ERRORS(sprintf('Adaptive stage %d',k),stages{k});
    F{k} = MAP_SHOW_SINGLE('fourleaf',stages{k}.xiMAP,d.angles, ...
        'stage',sprintf('FOURLEAF | Experiment 1 | Stage %d',k),k);
end
R.adaptiveStages = stages;
R.adaptiveFinal = stages{4};
out = struct('shape',S.kind,'experiment',1, ...
    'noiseCase',S.noiseCase,'replicate',S.rep, ...
    'route',S.route,'noiseSeed',S.noise.seed, ...
    'designInfo',S.designInfo,'fits',R,'figures',{F});
% No save, saveas, savefig, print, or exportgraphics calls.
fprintf('All figures remain open. Save individual figures manually.\n');
end
