%% Paper tables - from the result files of Simulated_groundtruth_evaluation and Real_datasets_evaluation
%
% Reads results/*.csv and prints and writes:
%   1. final_tables_    mean and SD over recordings of every metric per dataset (FINAL_*.csv), percentage of recordings that exceed each
%                       objective indicator (PREV_*.csv; ERP_distortion_dementia.csv for the participants with a distorted ERP), Spearman correlation of the pdBSI before/after denoising, ERP distortion in Dementia
%   2. final_p300_sim_  P300 of the healthy participants and the numbers quoted from the simulation (console)
% Run Simulated_groundtruth_evaluation and Real_datasets_evaluation (all datasets) first.

addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');
final_tables_();
final_p300_sim_();

function final_tables_()
% Final per-dataset summary tables (mean +- SD over recordings) and prevalence of the objective indicators.
rd = cfg_path('results');
cols = {'RMS_ratio','BPR_delta','BPR_theta','BPR_alpha','BPR_beta','Kurtosis','CovPreserve','ACF_corr','Entropy','GA50','GA100','GA150','pdBSI_delta','pdBSI_theta','pdBSI_alpha','pdBSI_beta','pdBSI_pre_delta','pdBSI_pre_theta','pdBSI_pre_alpha','pdBSI_pre_beta'};
% ---------------------------------------------------------------- Stroke (paper outputs, 718)
T = readtable(fullfile(rd, 'clin2_stroke_saved.csv')); T.Method = string(T.Method); T.Rec = string(T.Rec);
summ(T, cols, fullfile(rd, 'FINAL_stroke_paper.csv'), 'Stroke, saved denoised files');
prev(T, 'Stroke (718)', fullfile(rd, 'PREV_stroke.csv'));
% Spearman of pdBSI before/after
ms = unique(T.Method, 'stable')'; fprintf('\nStroke pdBSI per band: mean before / after, Spearman across recordings\n');
for m = ms(ms ~= "None")
    s = T(T.Method == m, :); fprintf('%-12s', m);
    for b = ["delta","theta","alpha","beta"]
        x = s.("pdBSI_pre_" + b); y = s.("pdBSI_" + b); ok = isfinite(x) & isfinite(y);
        fprintf(' %s: %.3f->%.3f rho=%.2f |', b, mean(x(ok)), mean(y(ok)), corr(x(ok), y(ok), 'Type', 'Spearman'));
    end
    fprintf('\n');
end
% ---------------------------------------------------------------- Stroke sensitivity (Original files, 718)
T2 = readtable(fullfile(rd, 'clin2_stroke.csv')); T2.Method = string(T2.Method);
summ(T2, cols, fullfile(rd, 'FINAL_stroke_sens.csv'), 'Stroke, defaults and sensitivity variants (Original files)');
% ---------------------------------------------------------------- Dementia (24)
T3 = readtable(fullfile(rd, 'clin2_dementia.csv')); T3.Method = string(T3.Method); T3.Rec = string(T3.Rec);
summ(T3, cols, fullfile(rd, 'FINAL_dementia.csv'), 'Dementia (24)');
% ---------------------------------------------------------------- AutomaticArtefactRemoval (13)
T4 = readtable(fullfile(rd, 'clin_auto.csv')); T4.Method = string(T4.Method); T4.Rec = string(T4.Rec);
for nc = [27 16 8 4]
    S = T4(T4.Nch == nc, :); summ(S, cols, fullfile(rd, sprintf('FINAL_auto_%dch.csv', nc)), sprintf('AutomaticArtefactRemoval, %d channels (13 participants)', nc));
end
prev(T4(T4.Nch == 27, :), 'Auto 27ch', fullfile(rd, 'PREV_auto27.csv')); prev(T4(T4.Nch == 4, :), 'Auto 4ch', fullfile(rd, 'PREV_auto4.csv'));
% ---------------------------------------------------------------- Dementia prevalence incl. ERP distortion
P = readtable(fullfile(rd, 'p300_subject_level.csv')); P.Method = string(P.Method); P.Rec = string(P.Rec);
prev(T3, 'Dementia (24)', fullfile(rd, 'PREV_dementia.csv'));
fprintf('\nDementia: participants with ERP waveform correlation < 0.8 (mean over 3 conditions x 2 regions), by method\n');
erp = table();
for m = unique(P.Method, 'stable')'
    if m == "None", continue; end
    s = P(P.Method == m, :); r = unique(s.Rec); v = arrayfun(@(i) mean(s.ERP_corr(s.Rec == r(i))), 1:numel(r));
    fprintf('%-10s %2d of %d participants (%.0f%%)\n', m, sum(v < 0.8), numel(r), 100*mean(v < 0.8));
    erp = [erp; table(m, numel(r), sum(v < 0.8), 100*mean(v < 0.8), 'VariableNames', {'Method', 'N', 'N_ERP_corr_below_0p8', 'pct_ERP_corr_below_0p8'})]; %#ok<AGROW>
end
writetable(erp, fullfile(rd, 'ERP_distortion_dementia.csv'));
end

function final_p300_sim_()
rd = cfg_path('results');
% ---------------- P300, controls only (9)
L = load(fullfile(rd, 'p300_erps2.mat')); ids = string(L.ids); grp = string(L.grp); ctrl = ids(grp == "C");
P = readtable(fullfile(rd, 'p300_subject_level.csv')); P.Method = string(P.Method); P.Rec = string(P.Rec); P = P(ismember(P.Rec, ctrl), :);
G = readtable(fullfile(rd, 'p300_group_level.csv')); G.Method = string(G.Method); G = G(string(G.Group) == "C", :);
fprintf('P300, controls (n=%d): group level = grand average of the controls; SNR and trial SD = median over controls x conditions x regions\n', numel(ctrl));
ms = unique(P.Method, 'stable')'; b = P(P.Method == "None", :);
for m = ms
    s = P(P.Method == m, :); g = G(G.Method == m, :);
    [~, ia, ib] = intersect(strcat(s.Rec, s.Condition, s.Region), strcat(b.Rec, b.Condition, b.Region));
    sdr = median(s.TrialSD_uV(ia) ./ b.TrialSD_uV(ib), 'omitnan'); snr = median(s.ERP_SNR, 'omitnan'); snr0 = median(b.ERP_SNR, 'omitnan');
    fprintf('%-10s waveform r=%.2f | grand-avg mean amp=%.2f uV [%.2f, %.2f] | peak latency=%.0f ms [%.0f, %.0f] | ERP SNR median=%.1f (None %.1f) | trial SD ratio=%.2f\n', m, mean(g.Waveform_corr_vs_None), mean(g.MeanAmp_uV), mean(g.MeanAmp_lo), mean(g.MeanAmp_hi), mean(g.PeakLat_ms), mean(g.PeakLat_lo), mean(g.PeakLat_hi), snr, snr0, sdr);
end
% ---------------- Simulation (final parameters)
A = readtable(fullfile(rd, 'SUMMARY_scale1.csv')); A.Method = string(A.Method);
V = readtable(fullfile(rd, 'sim_wt_coif4_variants.csv')); V.Method = string(V.Method);
fprintf('\nSIM, SNR (affected channels), mean (sd), by channel count\n');
for m = ["None","WT_repo","ASR_k5","ASR_k10","ASR_k20","ASR_k30","EMD_2IMF","EMD_adapt","ASR5_WT","EMD2_WT","WT_EMD2"]
    fprintf('%-10s', m); for c = [4 8 16 31], r = A(A.Method == m & A.Nch == c, :); fprintf('  %6.2f(%.2f)', r.SNR_aff_mean, r.SNR_aff_sd); end; fprintf('\n');
end
fprintf('\nSIM, RRMSE_t / RRMSE_f / CC at 31 and 4 ch\n');
for m = ["None","WT_repo","ASR_k5","ASR_k20","EMD_2IMF","EMD_adapt","ASR5_WT","EMD2_WT"]
    fprintf('%-10s', m); for c = [31 4], r = A(A.Method == m & A.Nch == c, :); fprintf('  [%d] RRt=%.2f RRf=%.2f CC=%.2f', c, r.RRMSE_t_aff_mean, r.RRMSE_f_aff_mean, r.CC_aff_mean); end; fprintf('\n');
end
fprintf('\nSIM, WT variants (coif4), SNR mean by channel count\n');
for m = unique(V.Method, 'stable')'
    fprintf('%-14s', m); for c = [4 8 16 31], s = V(V.Method == m & V.Nch == c, :); fprintf('  %6.2f', mean(s.SNR_aff)); end; fprintf('\n');
end
fprintf('\nSIM other variants (31ch / 4ch SNR)\n');
for m = ["WT_repo_bior44","WT_repo_db4","WT_repo_sym4","WT_den_hard","WT_den_soft","EMD_3IMF","EMD_4IMF","EMD_rm2","WT_EMD2"]
    r31 = A(A.Method == m & A.Nch == 31, :); r4 = A(A.Method == m & A.Nch == 4, :); fprintf('%-16s %6.2f / %6.2f\n', m, r31.SNR_aff_mean, r4.SNR_aff_mean);
end
B = readtable(fullfile(rd, 'SUMMARY_scales.csv')); B.Method = string(B.Method);
fprintf('\nSIM SNR by artefact scale (x0.25 x0.5 x1 x2 x4)\n');
for c = [31 4], for m = ["None","WT_repo","ASR_k5","ASR_k20"], fprintf('%2dch %-8s', c, m); for sc = [0.25 0.5 1 2 4], r = B(B.Method == m & B.Nch == c & B.Scale == sc, :); fprintf(' %6.2f', r.SNR_aff); end; fprintf('\n'); end, end
end

function summ(T, cols, fname, label)
ms = unique(T.Method, 'stable')'; rows = {};
for m = ms
    s = T(T.Method == m, :); r = {char(m), height(s)};
    for k = 1:numel(cols), v = s.(cols{k}); r{end+1} = mean(v, 'omitnan'); r{end+1} = std(v, 'omitnan'); end %#ok<AGROW>
    rows(end+1,:) = r; %#ok<AGROW>
end
vn = {'Method','N'}; for k = 1:numel(cols), vn = [vn, {[cols{k} '_mean'], [cols{k} '_sd']}]; end %#ok<AGROW>
A = cell2table(rows, 'VariableNames', vn); writetable(A, fname);
fprintf('\n=== %s (N=%d recordings)\n', label, A.N(1));
disp(A(:, {'Method','RMS_ratio_mean','BPR_delta_mean','BPR_theta_mean','BPR_alpha_mean','Kurtosis_mean','CovPreserve_mean','ACF_corr_mean','GA50_mean','GA100_mean'}));
end

function prev(T, label, fname)
ms = unique(T.Method, 'stable')'; rows = {};
for m = ms
    s = T(T.Method == m, :); n = height(s);
    rows(end+1,:) = {char(m), n, 100*mean(s.GA100 < 95), 100*mean(s.RMS_ratio < 0.5), 100*mean(s.BPR_delta < 0.5)}; %#ok<AGROW>
end
A = cell2table(rows, 'VariableNames', {'Method','N','pct_GA100_below95','pct_RMS_below0p5','pct_delta_below50'}); writetable(A, fname);
fprintf('\n--- Prevalence of indicators, %s\n', label); disp(A);
end
