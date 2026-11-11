function [out, files] = RUN_AND_SAVE_PAPER(entryPoint,noiseCase,rep,outputDir)
%RUN_AND_SAVE_PAPER Run an UNMODIFIED experiment entry and save its results.
% The original *_EXP*.m functions remain the authoritative experiments.
% Plotting and figure handling remain exactly as in the original scripts.
%
% Usage:
%   [out,files] = RUN_AND_SAVE_PAPER('CIRCLE_EXP1');
%   [out,files] = RUN_AND_SAVE_PAPER('PEANUT_EXP3','noisy',1);
%   [out,files] = RUN_AND_SAVE_PAPER('KITE_EXP2','noisy',1,'my-results');
if nargin<2 || isempty(noiseCase), noiseCase='noisy';end
if nargin<3 || isempty(rep), rep=1;end
if nargin<4, outputDir=[];end
entryPoint=upper(char(entryPoint));
shapes={'CIRCLE','KITE','FOURLEAF','PEANUT'};
valid=false;
for j=1:numel(shapes)
    for e=1:3
        if strcmp(entryPoint,sprintf('%s_EXP%d',shapes{j},e))
            valid=true;
        end
    end
end
assert(valid,'Unknown experiment entry point: %s',entryPoint);
% No code path for loading old MAP estimates: reproduce from scratch.
out=feval(str2func(entryPoint),noiseCase,rep);
files=SAVE_PAPER_RESULT(out,outputDir);
end
