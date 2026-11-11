function summary = RUN_ALL_AND_SAVE_PAPER(noiseCase,rep,outputDir)
%RUN_ALL_AND_SAVE_PAPER Recompute all 12 paper experiments and save results.
%   RUN_ALL_AND_SAVE_PAPER
%   RUN_ALL_AND_SAVE_PAPER('noisy',1)
%   S = RUN_ALL_AND_SAVE_PAPER('noisy',1,'my_results');
%
% This is an ADD-ON runner: it does not change any original experiment,
% plotting routine, inverse solver, or fixed-noise file. All experiments
% run from scratch. Figure windows stay open as in the original code.
% The numerical results are saved separately as MAT and CSV files, plus
% a human-readable TXT summary. Existing results are never overwritten.
if nargin < 1 || isempty(noiseCase), noiseCase = 'noisy'; end
if nargin < 2 || isempty(rep), rep = 1; end
root = fileparts(mfilename('fullpath'));
if nargin < 3 || isempty(outputDir)
    outputDir = fullfile(root,'results');
end
noiseCase = lower(char(noiseCase));
outputDir = char(outputDir);
assert(any(strcmp(noiseCase,{'noisy','zero'})), ...
    'noiseCase must be ''noisy'' or ''zero''.');
assert(isscalar(rep) && rep >= 1 && rep == floor(rep), ...
    'rep must be a positive integer.');
CHECK_PROJECT();
if exist(outputDir,'dir') ~= 7
    [ok,msg] = mkdir(outputDir);
    assert(ok,'Cannot create results folder: %s',msg);
end

% Use a new name every time, just like SAVE_PAPER_RESULT.
base = sprintf('ALL_12_%s_rep%d',noiseCase,rep);
runIndex = 1;
while true
    txtFile = fullfile(outputDir,sprintf('%s_run%03d.txt',base,runIndex));
    if exist(txtFile,'file') ~= 2
        break;
    end
    runIndex = runIndex + 1;
end
fid = fopen(txtFile,'wt');
assert(fid > 0,'Cannot write summary text file: %s',txtFile);
fprintf(fid,'Moving-sensor inverse source problem: all 12 MATLAB experiments\n');
fprintf(fid,'Noise case: %s, replicate: %d\n',noiseCase,rep);
fprintf(fid,'Generated on: %s\n',datestr(now));
fprintf(fid,'All MAP estimates are recomputed, not loaded from saved data.\n');
fprintf(fid,'Numerical MAT and CSV files are listed below.\n\n');

shapes = {'CIRCLE','KITE','FOURLEAF','PEANUT'};
records = repmat(struct('experiment','','mat','','csv','','completed',false),12,1);
count = 0;
try
    for s = 1:numel(shapes)
        for e = 1:3
            count = count + 1;
            entry = sprintf('%s_EXP%d',shapes{s},e);
            fprintf('\n======== %d/12: %s ========\n',count,entry);
            fprintf(fid,'======== %d/12: %s ========\n',count,entry);
            % Reuse the existing single-experiment runner unchanged.
            [out,files] = RUN_AND_SAVE_PAPER(entry,noiseCase,rep,outputDir);
            records(count).experiment = entry;
            records(count).mat = files.mat;
            records(count).csv = files.csv;
            records(count).completed = true;
            fprintf(fid,'MAT: %s\nCSV: %s\n',files.mat,files.csv);
            grid=round(out.route.mainWindows/(1/400));
            fprintf(fid,'Observation windows (grid indices): %s\n',mat2str(grid));
            localWriteFits(fid,out.fits);
            fprintf(fid,'\n');
            fflushIfAvailable(fid);
            clear out;
        end
    end
    fprintf(fid,'All 12 experiments completed successfully.\n');
    fclose(fid);
catch ME
    fprintf(fid,'\nFAILED at experiment %d (%s): %s\n',count,entry,ME.message);
    fclose(fid);
    fprintf(2,'Stopped on %s. Completed cases are already saved.\n',entry);
    fprintf(2,'Partial TXT summary: %s\n',txtFile);
    rethrow(ME);
end
summary = struct('results',records,'textFile',txtFile, ...
    'outputDir',outputDir,'noiseCase',noiseCase,'replicate',rep);
fprintf('\nAll 12 experiments finished.\nTXT summary: %s\n',txtFile);
fprintf('Numerical MAT/CSV files: %s\n',outputDir);
fprintf('MATLAB figure windows remain open; figures are not auto-exported.\n');
end

function localWriteFits(fid,fits)
fields = fieldnames(fits);
for j = 1:numel(fields)
    policy = fields{j};
    if strcmp(policy,'adaptiveStages')
        for k = 1:numel(fits.adaptiveStages)
            localWriteOne(fid,sprintf('adaptive_stage%d',k), ...
                fits.adaptiveStages{k});
        end
    elseif strcmp(policy,'adaptiveFinal') && isfield(fits,'adaptiveStages')
        % Already included as the last stage; avoid duplicate rows.
        continue;
    else
        localWriteOne(fid,policy,fits.(policy));
    end
end
end

function localWriteOne(fid,policy,fit)
x = fit.xiMAP(:);
xiString = strtrim(sprintf('%.12g ',x));
fprintf(fid,['  %-19s nMain=%3d  xi=[%s]  E_xi=%.12g  ', ...
    'E_D=%.12g  IoU=%.12g  MAPobjective=%.12g\n'], ...
    policy,fit.nMainObservations,xiString, ...
    fit.relativeParameterError,fit.relativeSymmetricDifference, ...
    fit.IoU,fit.objective);
end

function fflushIfAvailable(fid)
% Flush the summary periodically where supported; not required by MATLAB.
if exist('fflush','builtin') == 5 || exist('fflush','file') == 2
    fflush(fid);
end
end
