%% Computation times: mean and SD of the elapsed time, CPU time and memory of every method
%
% Real datasets: reads the files denoising_performance_*.csv written by Dataset_Stroke_new.m, Dataset_Dementia.m and
%   Dataset_automaticArtifactRemoval.m (one row per recording and method). Each recording is counted once.
% Semi-simulated data: reads results/sim_infos.mat written by Simulated_groundtruth_evaluation.m. WT and EMD are timed on the whole
%   31-channel recording (WT_repo, EMD, EMD2_WT = WT after the EMD); ASR is timed for each number of channels (ASR_k<k>_ch<n>).
%   Only the tasks computed in the last run are in the file.
% The times are recorded with the MATLAB profiler active (real datasets) and are indicative, not controlled benchmarks.
% Outputs (results/): computation_times_real.csv, computation_times_simulated.csv

addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');
rd = cfg_path('results');

%% Real datasets
roots = {'Stroke', cfg_path('stroke'); ...
         'Dementia', fullfile(cfg_path('dementia'), 'EEG_denoised'); ...
         'AutomaticArtefactRemoval', fullfile(cfg_path('auto'), 'Results')};
labels = containers.Map({'W', 'ASR', 'EMD', 'EMDW', 'ASRW'}, {'WT', 'ASR', 'EMD', 'EMD->WT', 'ASR->WT'});
rows = {};
for r = 1:size(roots, 1)
    files = dir(fullfile(roots{r, 2}, '**', 'denoising_performance_*.csv'));
    if isempty(files), fprintf('%s: no denoising_performance file in %s\n', roots{r, 1}, roots{r, 2}); continue; end
    T = table();
    for i = 1:numel(files)
        t = readtable(fullfile(files(i).folder, files(i).name));
        t.Subject = string(t.Subject); t.Method = string(t.Method);
        T = [T; t(:, {'Subject', 'Method', 'ElapsedTime', 'CPUTime', 'MemoryUsed'})]; %#ok<AGROW>
    end
    [~, keep] = unique(T.Subject + "|" + T.Method, 'stable'); T = T(keep, :);   % files of the same recording can repeat rows
    tok = regexp(T.Method, '^(EEG_\d+)_(.*)$', 'tokens', 'once');
    config = repmat("", height(T), 1); method = T.Method;
    for i = 1:height(T)
        if ~isempty(tok{i}), config(i) = tok{i}{1}; method(i) = tok{i}{2}; end
    end
    for g = unique(config)'
        for m = ["W", "ASR", "EMD", "EMDW", "ASRW"]
            s = T(config == g & method == m, :);
            if isempty(s), continue; end
            rows(end+1, :) = {roots{r, 1}, char(g), labels(char(m)), height(s), mean(s.ElapsedTime), std(s.ElapsedTime), ...
                mean(s.CPUTime), std(s.CPUTime), mean(s.MemoryUsed), std(s.MemoryUsed)}; %#ok<AGROW>
        end
    end
end
if ~isempty(rows)
    R = cell2table(rows, 'VariableNames', {'Dataset', 'Config', 'Method', 'N', 'Elapsed_s_mean', 'Elapsed_s_sd', 'CPU_s_mean', 'CPU_s_sd', 'Memory_MB_mean', 'Memory_MB_sd'});
    writetable(R, fullfile(rd, 'computation_times_real.csv')); disp(R);
end

%% Semi-simulated data
f = fullfile(rd, 'sim_infos.mat');
if isfile(f)
    S = load(f, 'tms'); tms = S.tms(~cellfun(@isempty, S.tms));
    names = cellfun(@fieldnames, tms, 'UniformOutput', false); names = unique(vertcat(names{:}));
    rows = {};
    for k = 1:numel(names)
        v = cellfun(@(t) getfield_or_nan(t, names{k}), tms);
        v = v(~isnan(v));
        rows(end+1, :) = {names{k}, numel(v), mean(v), std(v)}; %#ok<AGROW>
    end
    R = cell2table(rows, 'VariableNames', {'Step', 'N', 'Elapsed_s_mean', 'Elapsed_s_sd'});
    writetable(R, fullfile(rd, 'computation_times_simulated.csv')); disp(R);
else
    fprintf('Semi-simulated: %s not found (run Simulated_groundtruth_evaluation first)\n', f);
end

function v = getfield_or_nan(t, name)
    if isfield(t, name), v = t.(name); else, v = NaN; end
end
