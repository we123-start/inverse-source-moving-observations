function out = CIRCLE_EXP2(noiseCase,rep)
%CIRCLE_EXP2 Independent Experiment 2 for the circle source.
% Original Q2 polar FEM + backward Euler + constrained SQP MAP.
% Usage: CIRCLE_EXP2                     (original fixed noisy realization)
%        CIRCLE_EXP2('zero',1)           (zero-noise experiment)
% Each call recomputes its OWN measurements/solutions from scratch.
% All figures stay OPEN; images and MAT results are NEVER auto-saved.
if nargin < 1 || isempty(noiseCase), noiseCase = 'noisy'; end
if nargin < 2 || isempty(rep), rep = 1; end
S = MAP_PREPARE_CASE('circle',noiseCase,rep);
cfg = S.cfg;
w0 = S.w0;
R = struct();
F = cell(1,4);
% Identical 320 main observation times, fixed noise vector and MAP start.
% The final adaptive solution is recomputed here, NOT loaded from Exp 1.
policies = {'fixed_initial','fixed_low_A','fixed_low_B'};
R.adaptiveFinal = MAP_SOLVE_MEMORY(S.invModel,'circle', ...
    S.data.adaptive,w0,cfg);
MAP_PRINT_ERRORS('Adaptive final',R.adaptiveFinal);
F{1} = MAP_SHOW_SINGLE('circle',R.adaptiveFinal.xiMAP, ...
    S.data.adaptive.angles,'adaptive_final_compare', ...
    'CIRCLE | Experiment 2 | Adaptive reference',1);
for j = 1:3
    name = policies{j};
    R.(name) = MAP_SOLVE_MEMORY(S.invModel,'circle', ...
        S.data.(name),w0,cfg);
    MAP_PRINT_ERRORS(name,R.(name));
    F{j+1} = MAP_SHOW_SINGLE('circle',R.(name).xiMAP, ...
        S.data.(name).angles,'policy', ...
        sprintf('CIRCLE | Experiment 2 | %s',name),j+1);
end
out = struct('shape',S.kind,'experiment',2, ...
    'noiseCase',S.noiseCase,'replicate',S.rep, ...
    'route',S.route,'noiseSeed',S.noise.seed, ...
    'designInfo',S.designInfo,'fits',R,'figures',{F});
% No save, saveas, savefig, print, or exportgraphics calls.
fprintf('All figures remain open. Save individual figures manually.\n');
end
