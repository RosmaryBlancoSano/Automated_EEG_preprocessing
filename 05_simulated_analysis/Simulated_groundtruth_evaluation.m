%% Semi-simulated data scored against the CLEAN signal (ground truth), with sensitivity analysis
%
% Data: generated here (fake_EEG.m + EOG/EMG epochs of EEGdenoiseNet, linear mixing); the clean signal is kept.
% Methods: W, EMD, ASR and the hybrids, default settings and parameter variants (see 05_simulated_analysis/groundtruth_task.m).
% Metrics: SNR, RRMSE (time, frequency) and correlation against the clean signal (07_metrics/metrics_groundtruth.m).
% Outputs (results/): sim_s*_x*.csv (one file per subject and artefact amplitude), sim_wt_coif4_variants.csv,
%          ALL_sim_rows.csv, SUMMARY_scale1.csv, SUMMARY_scales.csv, SUMMARY_friedman.csv.
%
% Before running: download EEGdenoiseNet (EOG_all_epochs.mat, EMG_all_epochs.mat) and set its folder in cfg_path.m.
% Optional variables set before running: n_subjects (default 10), n_workers (default 4; 0 = no parallel pool), amplitudes (default [0.25 0.5 2 4]).
% (Reproduce_simulated.m runs this script and then builds the tables and the figure.)

addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup'); outdir = cfg_path('results');
if ~exist('n_subjects', 'var'), n_subjects = 10; end
if ~exist('n_workers', 'var'), n_workers = 4; end
if ~exist('amplitudes', 'var'), amplitudes = [0.25 0.5 2 4]; end
fs = 256;

%% 1. Run: n_subjects subjects (seeds 1..n) at artefact amplitude x1 with all variants, and at x0.25, x0.5, x2, x4 with the main methods
tasks = [];  % [seed amplitude full(1)/main(0)]
for s = 1:n_subjects, tasks(end+1,:) = [s 1 1]; end %#ok<SAGROW>
for sc = amplitudes, for s = 1:n_subjects, tasks(end+1,:) = [s sc 0]; end, end
start_pool(n_workers);
n = size(tasks, 1); infos = cell(n, 1); tms = cell(n, 1); t0 = tic;
parfor (i = 1:n, n_workers)
    seed = tasks(i,1); sc = tasks(i,2); lay = 'main'; if tasks(i,3) == 1, lay = 'full'; end
    f = fullfile(outdir, sprintf('sim_s%02d_x%g_%s.csv', seed, sc, lay));
    if ~isfile(f)
        try
            [T, info, tm] = groundtruth_task(seed, sc, lay);
            writetable(T, f); infos{i} = info; tms{i} = tm;
        catch ME
            fid = fopen(fullfile(outdir, sprintf('ERR_s%02d_x%g_%s.txt', seed, sc, lay)), 'w'); fprintf(fid, '%s\n', getReport(ME)); fclose(fid);
        end
    end
end
save(fullfile(outdir, 'sim_infos.mat'), 'infos', 'tasks', 'tms');
fprintf('Run finished in %.0f s\n', toc(t0));

%% 2. WT variants with the Coiflet 4 wavelet of the paper (level 5 and 9, soft rule, universal threshold), amplitude x1
res = cell(n_subjects, 1);
parfor (s = 1:n_subjects, n_workers)
    [clean, contam, info] = sim_make(s, 1); Cf = common_filter(clean, fs); Xf = common_filter(contam, fs);
    Channels = info.channels; ix = @(names) find(ismember(Channels, names));
    ex16 = {'AFz','Fz','FT7','FC3','FCz','FC4','FT8','Cz','TP7','CP3','CPz','CP4','TP8','Pz','OZ'};
    c16 = Channels(~ismember(Channels, ex16)); c8 = c16(~ismember(c16, {'F7','F3','F4','F8','P7','P3','P4','P8'})); c4 = c8(~ismember(c8, {'C3','C4','O1','O2'}));
    cfgs = struct('n', {31, 16, 8, 4}, 'rows', {1:31, ix(c16), ix(c8), ix(c4)});
    V = {'WT_coif4_L5', 'coif4', 5, 'Hard', 'Bayes'; 'WT_coif4_L9', 'coif4', 9, 'Hard', 'Bayes'; 'WT_coif4_soft', 'coif4', 7, 'Soft', 'Bayes'; 'WT_coif4_univ', 'coif4', 7, 'Hard', 'UniversalThreshold'};
    rows = {};
    for v = 1:size(V, 1)
        Y = wt_apply(Xf, V{v,2}, V{v,3}, V{v,4}, 'resid', V{v,5});
        for c = 1:numel(cfgs)
            M = metrics_groundtruth(Cf, Y, Xf, fs, cfgs(c).rows); a = ismember(cfgs(c).rows, info.affected);
            rows(end+1,:) = {s, cfgs(c).n, V{v,1}, mean(M.SNR(a)), mean(M.RRMSE_t(a)), mean(M.RRMSE_f(a)), mean(M.CC(a))}; %#ok<PFBNS>
        end
    end
    res{s} = cell2table(rows, 'VariableNames', {'Subject','Nch','Method','SNR_aff','RRMSE_t_aff','RRMSE_f_aff','CC_aff'});
end
T = vertcat(res{:}); writetable(T, fullfile(outdir, 'sim_wt_coif4_variants.csv'));

%% 3. Aggregate: summaries over subjects, Friedman tests, ranking at 31 and 4 channels
rd = outdir;
fs_ = dir(fullfile(rd, 'sim_s*_x*.csv')); T = table();
for i = 1:numel(fs_), T = [T; readtable(fullfile(rd, fs_(i).name))]; end %#ok<AGROW>
T.Method = string(T.Method);
writetable(T, fullfile(rd, 'ALL_sim_rows.csv'));
S1 = T(T.Scale == 1, :);
mets = {'SNR','RRMSE_t','RRMSE_f','CC','SNR_aff','RRMSE_t_aff','RRMSE_f_aff','CC_aff','SNR_repo','RRMSE_t_repo','CC_repo'};
meths = unique(S1.Method, 'stable'); nchs = [4 8 16 31];
rows = {};
for m = 1:numel(meths)
    for c = nchs
        s = S1(S1.Method == meths(m) & S1.Nch == c, :);
        if isempty(s), continue; end
        r = {char(meths(m)), c, height(s)};
        for k = 1:numel(mets), v = s.(mets{k}); r{end+1} = mean(v, 'omitnan'); r{end+1} = std(v, 'omitnan'); end %#ok<AGROW>
        rows(end+1,:) = r; %#ok<AGROW>
    end
end
vn = {'Method','Nch','N'};
for k = 1:numel(mets), vn = [vn, {[mets{k} '_mean'], [mets{k} '_sd']}]; end %#ok<AGROW>
A = cell2table(rows, 'VariableNames', vn); writetable(A, fullfile(rd, 'SUMMARY_scale1.csv'));

% artefact amplitude (signal-to-artefact ratio) for the main methods
mm = ["None","WT_repo","EMD_2IMF","EMD_adapt","ASR_k5","ASR_k20","ASR5_WT","EMD2_WT"];
scs = unique(T.Scale)'; rows = {};
for m = mm
    for sc = scs
        for c = nchs
            s = T(T.Method == m & T.Scale == sc & T.Nch == c, :);
            if isempty(s), continue; end
            rows(end+1,:) = {char(m), sc, c, height(s), mean(s.SNR_aff), std(s.SNR_aff), mean(s.RRMSE_t_aff), std(s.RRMSE_t_aff), mean(s.CC_aff), std(s.CC_aff)}; %#ok<AGROW>
        end
    end
end
B = cell2table(rows, 'VariableNames', {'Method','Scale','Nch','N','SNR_aff','SNR_sd','RRMSE_t_aff','RRMSE_t_sd','CC_aff','CC_sd'});
writetable(B, fullfile(rd, 'SUMMARY_scales.csv'));

% Friedman tests across the default methods at amplitude x1 (subjects = blocks)
dm = ["WT_repo","EMD_2IMF","ASR_k5","EMD2_WT","ASR5_WT"]; ft = {};
for c = nchs
    for k = {'SNR_aff','RRMSE_t_aff','CC_aff'}
        M = nan(n_subjects, numel(dm));
        for j = 1:numel(dm), s = S1(S1.Method == dm(j) & S1.Nch == c, :); [~, o] = sort(s.Subject); M(1:height(s), j) = s.(k{1})(o); end
        M = M(all(isfinite(M), 2), :);
        if size(M, 1) >= 3, p = friedman(M, 1, 'off'); else, p = NaN; end
        [~, best] = max(mean(M, 1) * (1 - 2*strcmp(k{1}, 'RRMSE_t_aff')));
        ft(end+1,:) = {c, k{1}, size(M, 1), p, char(dm(best))}; %#ok<AGROW>
    end
end
F = cell2table(ft, 'VariableNames', {'Nch','Metric','N','Friedman_p','Best_default_method'}); writetable(F, fullfile(rd, 'SUMMARY_friedman.csv'));

for c = [31 4]
    s = A(A.Nch == c, :); [~, o] = sort(s.SNR_aff_mean, 'descend');
    fprintf('\nRanking by SNR on the contaminated channels, %d channels\n', c);
    disp(s(o, {'Method','N','SNR_aff_mean','SNR_aff_sd','RRMSE_t_aff_mean','RRMSE_f_aff_mean','CC_aff_mean'}));
end
disp(F);
if isfile(fullfile(rd, 'sim_infos.mat'))
    L = load(fullfile(rd, 'sim_infos.mat')); ok = ~cellfun(@isempty, L.infos); ii = L.infos(ok); tk = L.tasks(ok, :);
    sar = []; for i = 1:numel(ii), if tk(i,2) == 1, sar(end+1) = median(ii{i}.SAR_dB); end, end %#ok<AGROW>
    if ~isempty(sar), fprintf('Median signal-to-artefact ratio (x1) across subjects: %.1f dB (range %.1f to %.1f)\n', median(sar), min(sar), max(sar)); end
end
