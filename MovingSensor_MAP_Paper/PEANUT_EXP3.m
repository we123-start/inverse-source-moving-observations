function out = PEANUT_EXP3(noiseCase,rep)
%PEANUT_EXP3 Independent Experiment 3 for the peanut source.
% Original Q2 polar FEM + backward Euler + constrained SQP MAP.
% Usage: PEANUT_EXP3                     (original fixed noisy realization)
%        PEANUT_EXP3('zero',1)           (zero-noise experiment)
% Each call recomputes its OWN measurements/solutions from scratch.
% All figures stay OPEN; images and MAT results are NEVER auto-saved.
if nargin < 1 || isempty(noiseCase), noiseCase = 'noisy'; end
if nargin < 2 || isempty(rep), rep = 1; end
S = MAP_PREPARE_CASE('peanut',noiseCase,rep);
cfg = S.cfg;
w0 = S.w0;
R = struct();
F = cell(1,4);
% Two-site, predetermined nonadaptive moving controls.
% Low-flux sites are truth-informed synthetic stress tests, not an
% implementable blind site-selection method.
% Recompute adaptive final locally: NO dependency on Experiment 1.
policies = {'local_A','local_B','local_C'};
R.adaptiveFinal = MAP_SOLVE_MEMORY(S.invModel,'peanut', ...
    S.data.adaptive,w0,cfg);
MAP_PRINT_ERRORS('Adaptive final',R.adaptiveFinal);
F{1} = MAP_SHOW_SINGLE('peanut',R.adaptiveFinal.xiMAP, ...
    S.data.adaptive.angles,'adaptive_final_compare', ...
    'PEANUT | Experiment 3 | Adaptive reference',1);
for j = 1:3
    name = policies{j};
    R.(name) = MAP_SOLVE_MEMORY(S.invModel,'peanut', ...
        S.data.(name),w0,cfg);
    MAP_PRINT_ERRORS(name,R.(name));
    F{j+1} = MAP_SHOW_SINGLE('peanut',R.(name).xiMAP, ...
        S.data.(name).angles,'comparison_moving', ...
        sprintf('PEANUT | Experiment 3 | %s',name),j+1);
end
out = struct('shape',S.kind,'experiment',3, ...
    'noiseCase',S.noiseCase,'replicate',S.rep, ...
    'route',S.route,'noiseSeed',S.noise.seed, ...
    'designInfo',S.designInfo,'fits',R,'figures',{F});
% No save, saveas, savefig, print, or exportgraphics calls.
fprintf('All figures remain open. Save individual figures manually.\n');
end
