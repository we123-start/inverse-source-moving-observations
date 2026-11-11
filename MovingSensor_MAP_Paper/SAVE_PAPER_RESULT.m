function files = SAVE_PAPER_RESULT(out, outputDir)
%SAVE_PAPER_RESULT Explicit, non-destructive saving of computed MAP results.
% This is an OPTIONAL utility. It does not change the original experiment
% scripts, numerical solver, plot style, figures, or default save behavior.
%
% Examples:
%   out = CIRCLE_EXP1;
%   files = SAVE_PAPER_RESULT(out);
%   files = SAVE_PAPER_RESULT(out,'D:\\my_results');
%
% The .mat file stores all original numerical output except live Figure
% handles, which cannot be used after MATLAB closes. The .csv file contains
% the parameter estimates and errors needed for the manuscript tables.
% Existing files are never overwritten: a fresh _runNNN suffix is added.
assert(isstruct(out) && all(isfield(out, ...
    {'shape','experiment','noiseCase','replicate','fits'})), ...
    'Expected a result returned by one of the twelve *_EXP* functions.');
if nargin < 2 || isempty(outputDir)
    outputDir = fullfile(fileparts(mfilename('fullpath')),'results');
end
outputDir = char(outputDir);
if exist(outputDir,'dir') ~= 7
    [ok,msg] = mkdir(outputDir);
    assert(ok,'Cannot create result directory %s: %s',outputDir,msg);
end
shape = upper(char(out.shape));
noiseCase = lower(char(out.noiseCase));
expNumber = out.experiment;
rep = out.replicate;
assert(any(strcmp(shape,{'CIRCLE','KITE','FOURLEAF','PEANUT'})));
assert(isscalar(expNumber) && any(expNumber==1:3));
assert(isscalar(rep) && rep>=1 && rep==floor(rep));
base = sprintf('%s_EXP%d_%s_rep%d',shape,expNumber,noiseCase,rep);
i = 1;
while true
    stem = sprintf('%s_run%03d',base,i);
    matFile = fullfile(outputDir,[stem,'.mat']);
    csvFile = fullfile(outputDir,[stem,'.csv']);
    if exist(matFile,'file')~=2 && exist(csvFile,'file')~=2
        break;
    end
    i = i+1;
end
% Never serialize Figure handles or touch currently displayed figures.
numericResult = out;
if isfield(numericResult,'figures')
    numericResult = rmfield(numericResult,'figures');
end
provenance = struct('experiment',sprintf('%s_EXP%d',shape,expNumber), ...
    'noiseCase',noiseCase,'replicate',rep, ...
    'matlabVersion',version,'savedOnLocalTime',datestr(now,30), ...
    'note','Explicit save; no cached inverse results are loaded by any experiment');
save(matFile,'numericResult','provenance','-v7');

fid = fopen(csvFile,'wt');
assert(fid>0,'Cannot open CSV for writing: %s',csvFile);
try
    fprintf(fid,['shape,experiment,policy,stage,nMainObservations,', ...
        'xi1,xi2,xi3,xi4,xi5,E_xi,E_D,IoU,MAP_objective,', ...
        'firstPassExitFlag,refinementExitFlag,refinementAccepted,', ...
        'finalCandidateExitFlag\n']);
    fields = fieldnames(out.fits);
    for j = 1:numel(fields)
        name = fields{j};
        if strcmp(name,'adaptiveStages')
            arr = out.fits.(name);
            for k = 1:numel(arr)
                localWriteFit(fid,out,sprintf('adaptive_stage%d',k),k,arr{k});
            end
        elseif strcmp(name,'adaptiveFinal') && isfield(out.fits,'adaptiveStages')
            % The final result already appears as the last adaptive stage.
            continue;
        else
            fit = out.fits.(name);
            localWriteFit(fid,out,name,4,fit);
        end
    end
catch ME
    fclose(fid);
    rethrow(ME);
end
fclose(fid);
files = struct('mat',matFile,'csv',csvFile);
fprintf('Saved numerical data (not figures):\n  %s\n  %s\n',matFile,csvFile);
end

function localWriteFit(fid,out,policy,stage,fit)
assert(all(isfield(fit,{'xiMAP','relativeParameterError', ...
    'relativeSymmetricDifference','IoU','objective','nMainObservations'})));
xi = nan(5,1);
x = fit.xiMAP(:);
assert(numel(x)<=5 && all(isfinite(x)));
xi(1:numel(x))=x;
firstExit=NaN;
if isfield(fit,'trials') && isstruct(fit.trials) && ...
        isfield(fit.trials,'exitflag') && ~isempty(fit.trials)
    firstExit=fit.trials(1).exitflag;
end
refineExit=NaN;refineAccepted=0;finalExit=firstExit;
if isfield(fit,'trials') && isstruct(fit.trials) && ~isempty(fit.trials)
    T=fit.trials(1);
    if isfield(T,'refinementExitFlag'),refineExit=T.refinementExitFlag;end
    if isfield(T,'refinementAccepted'),refineAccepted=double(T.refinementAccepted);end
    if isfield(T,'finalCandidateExitFlag'),finalExit=T.finalCandidateExitFlag;end
end
fprintf(fid,'%s,%d,%s,%d,%d,',out.shape,out.experiment, ...
    policy,stage,fit.nMainObservations);
fprintf(fid,'%.15g,%.15g,%.15g,%.15g,%.15g,',xi);
fprintf(fid,'%.15g,%.15g,%.15g,%.15g,%.15g,%.15g,%d,%.15g\n', ...
    fit.relativeParameterError,fit.relativeSymmetricDifference, ...
    fit.IoU,fit.objective,firstExit,refineExit,refineAccepted,finalExit);
end
