%% Stroke biomarkers before and after denoising (718 recordings)
%
% For the signal before denoising (sub-XXXX_OriginalInterpolated.set) and for the denoised files saved by Dataset_Stroke_new.m (W, ASR, EMD, EMD_W, ASR_W):
% relative power in the delta, theta, alpha and beta bands, delta/alpha and theta/alpha ratios of the relative powers, and BSI per band,
% on non-overlapping 10 s windows averaged over windows (07_metrics/stroke_features.m). The Spearman correlation across the recordings
% between the value before and after denoising quantifies whether the ranking of the patients is preserved.
% The pdBSI rows of the biomarker table come from Paper_tables.m.
% Optional: n_workers (default 4).
% Outputs (results/): stroke_features.csv (one row per recording and method), stroke_biomarkers_spearman.csv (one row per method)

addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');
if ~exist('n_workers', 'var'), n_workers = 4; end
root = cfg_path('stroke'); outdir = cfg_path('results'); recdir = fullfile(outdir, 'rec_stroke_features'); if ~isfolder(recdir), mkdir(recdir); end
ids = readlines(cfg_path('ids718')); ids = ids(strlength(ids) > 0); subs = cellfun(@(s) ['sub-' s(end-3:end)], cellstr(ids), 'UniformOutput', false);
start_pool(n_workers);
n = numel(subs); t0 = tic;
parfor (i = 1:n, n_workers)
    sub = subs{i}; f = fullfile(recdir, [sub '.csv']);
    if ~isfile(f)
        try
            r = fullfile(root, sub);
            tags = {'Original','W','ASR','EMD','EMD_W','ASR_W'}; names = {'None','WT_paper','ASR_paper','EMD_paper','EMD_W_paper','ASR_W_paper'};
            [X, fs] = read_eeglab_set(fullfile(r, [sub '_OriginalInterpolated.set'])); rows = cell(numel(tags), 1);
            for k = 1:numel(tags)
                if k == 1, Y = X; else
                    Y = read_eeglab_set(fullfile(r, [sub '_' tags{k} '_denoised.set']));
                    if strcmp(tags{k}, 'ASR'), [xc, lags] = xcorr(X(1,:), Y(1,:), 'coeff'); [~, im] = max(xc); Y = circshift(Y, lags(im), 2); end
                end
                F = stroke_features(Y, fs); T = struct2table(F); T.Rec = string(sub); T.Method = string(names{k}); rows{k} = T;
            end
            writetable(vertcat(rows{:}), f);
        catch ME
            fid = fopen(fullfile(recdir, [sub '.ERR.txt']), 'w'); fprintf(fid, '%s\n', getReport(ME)); fclose(fid);
        end
    end
end
cs = dir(fullfile(recdir, '*.csv')); res = cell(numel(cs), 1);
for i = 1:numel(cs), res{i} = readtable(fullfile(recdir, cs(i).name)); end
T = vertcat(res{:}); T.Rec = string(T.Rec); T.Method = string(T.Method); writetable(T, fullfile(outdir, 'stroke_features.csv'));
fprintf('stroke_features: %d of %d recordings, %.0f s\n', numel(cs), n, toc(t0));

%% Spearman correlation across the recordings, before vs after denoising
feat = {'RelDelta', 'RelTheta', 'RelAlpha', 'RelBeta', 'DeltaAlpha', 'ThetaAlpha'};
methods = {'WT_paper', 'ASR_paper', 'EMD_paper', 'EMD_W_paper', 'ASR_W_paper'};
B = T(T.Method == "None", :); rho = nan(numel(methods), numel(feat));
for m = 1:numel(methods)
    A = T(T.Method == methods{m}, :); [~, ia, ib] = intersect(A.Rec, B.Rec);
    for f = 1:numel(feat)
        x = B.(feat{f})(ib); y = A.(feat{f})(ia); ok = isfinite(x) & isfinite(y);
        rho(m, f) = corr(x(ok), y(ok), 'Type', 'Spearman');
    end
end
S = array2table(rho, 'VariableNames', feat); S.Method = string(methods)'; S = movevars(S, 'Method', 'Before', 1);
writetable(S, fullfile(outdir, 'stroke_biomarkers_spearman.csv')); disp(S);
