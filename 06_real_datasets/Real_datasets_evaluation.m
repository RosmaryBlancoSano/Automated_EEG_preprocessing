%% Real datasets: preservation of the signal, sensitivity analysis, biomarkers and P300
%
% There is no clean reference, so each method is compared with the unprocessed (filtered) signal:
%   GA at +-50/100/150 uV, band powers, RMS amplitude, kurtosis, spectral entropy, correlation structure between channels,
%   autocorrelation, pair-derived Brain Symmetry Index (pdBSI) per band  ->  07_metrics/metrics_preservation.m
%
% Choose the dataset before running (default 'dementia'):
%   dataset = 'dementia'      24 recordings (ID0013 excluded), default settings + parameter variants
%   dataset = 'stroke_saved'  718 Stroke recordings, evaluates the denoised files saved by Dataset_Stroke_new.m
%   dataset = 'stroke'        718 Stroke recordings, re-runs the methods and the parameter variants
%   dataset = 'AutomaticArtefactRemoval'   13 participants, whole recordings, 27/16/8/4 channels
%   dataset = 'p300'          Dementia: P300 after denoising (Neutral/Positive/Negative pictures; controls and patients)
% Optional: n_workers (default 4).
% Outputs (results/): clin2_dementia.csv, clin2_stroke_saved.csv, clin2_stroke.csv, clin_auto.csv,
%          p300_subject_level.csv, p300_group_level.csv, p300_group_difference.csv, p300_erps2.mat
% Stroke and Dementia recordings are patient data and are not in the repository: set their folders in cfg_path.m.

addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');
if ~exist('dataset', 'var'), dataset = 'dementia'; end
if ~exist('n_workers', 'var'), n_workers = 4; end
outdir = cfg_path('results');
start_pool(n_workers);
t0 = tic;

switch dataset
    case {'dementia', 'stroke'}
        recdir = fullfile(outdir, ['rec2_' dataset]); if ~isfolder(recdir), mkdir(recdir); end
        if strcmp(dataset, 'dementia')
            d = dir(fullfile(cfg_path('dementia'), 'ID*', '*preproce*.set')); d = d(~contains({d.folder}, 'ID0013'));   % 24 participants as in the paper
            files = fullfile({d.folder}, {d.name}); pairs = [1 4; 2 3]; chsel = 1:4;                                   % TP9-TP10, AF7-AF8
        else
            ids = readlines(cfg_path('ids718')); ids = ids(strlength(ids) > 0); files = cell(1, numel(ids));
            for i = 1:numel(ids), nm = char(ids(i)); n4 = nm(end-3:end); files{i} = fullfile(cfg_path('stroke'), ['sub-' n4], ['sub-' n4 '_OriginalInterpolated.set']); end
            pairs = [1 2; 3 4; 5 6; 7 8]; chsel = 1:8;                                                                  % FC3-FC4, CP3-CP4, FT7-FT8, TP7-TP8
        end
        n = numel(files);
        parfor (i = 1:n, n_workers)
            [~, nm] = fileparts(files{i}); f = fullfile(recdir, [nm '.csv']);
            if ~isfile(f)
                try
                    [X, fs] = read_eeglab_set(files{i}); X = X(chsel, :);
                    R = denoise_and_score(X, fs, pairs, nm, 'sens'); R.Duration_s = repmat(size(X, 2)/fs, height(R), 1);
                    writetable(R, f);
                catch ME
                    fid = fopen(fullfile(recdir, [nm '.ERR.txt']), 'w'); fprintf(fid, '%s\n', getReport(ME)); fclose(fid);
                end
            end
        end
        cs = dir(fullfile(recdir, '*.csv')); res = cell(numel(cs), 1);
        for i = 1:numel(cs), res{i} = readtable(fullfile(recdir, cs(i).name)); end
        T = vertcat(res{:}); writetable(T, fullfile(outdir, sprintf('clin2_%s.csv', dataset)));
        fprintf('%s: %d of %d recordings, %.0f s\n', dataset, numel(cs), n, toc(t0));

    case 'stroke_saved'
        root = cfg_path('stroke'); recdir = fullfile(outdir, 'rec_stroke_saved2'); if ~isfolder(recdir), mkdir(recdir); end
        ids = readlines(cfg_path('ids718')); ids = ids(strlength(ids) > 0); subs = cellfun(@(s) ['sub-' s(end-3:end)], cellstr(ids), 'UniformOutput', false);
        pairs = [1 2; 3 4; 5 6; 7 8]; n = numel(subs);
        parfor (i = 1:n, n_workers)
            sub = subs{i}; f = fullfile(recdir, [sub '.csv']);
            if ~isfile(f)
                try
                    r = fullfile(root, sub); [X, fs] = read_eeglab_set(fullfile(r, [sub '_OriginalInterpolated.set']));
                    outs = struct(); tags = {'W','ASR','EMD','EMD_W','ASR_W'}; names = {'WT_paper','ASR_paper','EMD_paper','EMD_W_paper','ASR_W_paper'};
                    for k = 1:numel(tags)
                        Y = read_eeglab_set(fullfile(r, [sub '_' tags{k} '_denoised.set']));
                        if strcmp(tags{k}, 'ASR'), [xc, lags] = xcorr(X(1,:), Y(1,:), 'coeff'); [~, im] = max(xc); Y = circshift(Y, lags(im), 2); end   % realign the saved ASR output with the signal
                        outs.(names{k}) = Y;
                    end
                    R = denoise_and_score(X, fs, pairs, sub, 'outs', outs); R.Duration_s = repmat(size(X, 2)/fs, height(R), 1); writetable(R, f);
                catch ME
                    fid = fopen(fullfile(recdir, [sub '.ERR.txt']), 'w'); fprintf(fid, '%s\n', getReport(ME)); fclose(fid);
                end
            end
        end
        cs = dir(fullfile(recdir, '*.csv')); res = cell(numel(cs), 1);
        for i = 1:numel(cs), res{i} = readtable(fullfile(recdir, cs(i).name)); end
        T = vertcat(res{:}); writetable(T, fullfile(outdir, 'clin2_stroke_saved.csv'));
        fprintf('stroke_saved: %d of %d, %.0f s\n', numel(cs), n, toc(t0));

    case 'AutomaticArtefactRemoval'
        base = cfg_path('auto'); recdir = fullfile(outdir, 'rec_auto'); if ~isfolder(recdir), mkdir(recdir); end
        d = dir(fullfile(base, 's*')); d = d([d.isdir]); subs = {d.name}; subs = subs(~cellfun(@isempty, regexp(subs, '^s\d+$')));
        n = numel(subs);
        parfor (i = 1:n, n_workers)
            sub = subs{i}; f = fullfile(recdir, [sub '.csv']);
            if isfile(f), continue; end
            try
                R = score_auto_subject(fullfile(base, sub, [sub '.set']), sub);
                writetable(R, f);
            catch ME
                fid = fopen(fullfile(recdir, [sub '.ERR.txt']), 'w'); fprintf(fid, '%s\n', getReport(ME)); fclose(fid);
            end
        end
        cs = dir(fullfile(recdir, '*.csv')); res = cell(numel(cs), 1);
        for i = 1:numel(cs), res{i} = readtable(fullfile(recdir, cs(i).name)); end
        T = vertcat(res{:}); writetable(T, fullfile(outdir, 'clin_auto.csv'));
        fprintf('auto: %d of %d participants, %.0f s\n', numel(cs), n, toc(t0));

    case 'p300'
        % Subject level (amplitude, latency, ERP SNR, trial variability, waveform correlation), group level (grand averages with
        % bootstrap confidence intervals) and group difference (controls C vs patients P). Markers are in the last column of the Muse csv.
        d = dir(fullfile(cfg_path('dementia'), 'ID*')); d = d([d.isdir]); d = d(~strcmp({d.name}, 'ID0013')); files = {}; ids = {}; grp = {};   % 24 participants as in the paper
        for i = 1:numel(d)
            g = dir(fullfile(d(i).folder, d(i).name, '*Muse_output.csv'));
            if ~isempty(g), [~, ix] = max([g.bytes]); files{end+1} = fullfile(g(ix).folder, g(ix).name); ids{end+1} = d(i).name; grp{end+1} = g(ix).name(1); end %#ok<SAGROW>
        end
        n = numel(files); Rs = cell(n, 1); Es = cell(n, 1);
        parfor (i = 1:n, n_workers)
            try, [Rs{i}, Es{i}] = p300_subject(files{i}, ids{i}); catch ME, fprintf('ERR %s %s\n', ids{i}, ME.message); end
        end
        T = vertcat(Rs{:}); writetable(T, fullfile(outdir, 'p300_subject_level.csv')); save(fullfile(outdir, 'p300_erps2.mat'), 'Es', 'ids', 'grp');
        ok = find(~cellfun(@isempty, Es)); t = Es{ok(1)}.t; w = t >= 300 & t <= 550;
        meths = setdiff(fieldnames(Es{ok(1)}), {'t'}, 'stable'); conds = {'Positive','Neutral','Negative'}; regs = {'Frontal','Temporal'}; G = {'C','P'};
        rng(1); nb = 2000;
        rows = {};   % group level
        for m = 1:numel(meths)
            for c = 1:3
                for r = 1:2
                    for gi = 1:3
                        if gi <= 2, sel = ok(strcmp(grp(ok), G{gi})); lab = G{gi}; else, sel = ok; lab = 'All'; end
                        A = cell2mat(cellfun(@(e) e.(meths{m}).(conds{c}).(regs{r}), Es(sel), 'UniformOutput', false));
                        A0 = cell2mat(cellfun(@(e) e.None.(conds{c}).(regs{r}), Es(sel), 'UniformOutput', false));
                        ga = mean(A, 1); [pk, ip] = max(ga(w)); tw = t(w); lat = tw(ip); ma = mean(ga(w)); ma0 = mean(mean(A0, 1) .* w) / mean(w);
                        bm = zeros(nb, 1); bl = zeros(nb, 1);
                        for b = 1:nb
                            ix = randi(numel(sel), numel(sel), 1); g2 = mean(A(ix, :), 1); [~, ip2] = max(g2(w)); bl(b) = tw(ip2); bm(b) = mean(g2(w));
                        end
                        rows(end+1,:) = {meths{m}, conds{c}, regs{r}, lab, numel(sel), ma, prctile(bm, 2.5), prctile(bm, 97.5), pk, lat, prctile(bl, 2.5), prctile(bl, 97.5), corr(ga(:), mean(A0, 1)'), ma/ma0}; %#ok<AGROW>
                    end
                end
            end
        end
        GL = cell2table(rows, 'VariableNames', {'Method','Condition','Region','Group','nSubjects','MeanAmp_uV','MeanAmp_lo','MeanAmp_hi','Peak_uV','PeakLat_ms','PeakLat_lo','PeakLat_hi','Waveform_corr_vs_None','MeanAmp_ratio_vs_None'});
        writetable(GL, fullfile(outdir, 'p300_group_level.csv'));
        rows = {}; Cs = ok(strcmp(grp(ok), 'C')); Ps = ok(strcmp(grp(ok), 'P'));   % group difference, mean amplitude 300-550 ms
        for m = 1:numel(meths)
            for c = 1:3
                for r = 1:2
                    ac = arrayfun(@(i) mean(Es{i}.(meths{m}).(conds{c}).(regs{r})(w)), Cs); ap = arrayfun(@(i) mean(Es{i}.(meths{m}).(conds{c}).(regs{r})(w)), Ps);
                    pv = ranksum(ac, ap); pooled = sqrt(((numel(ac)-1)*var(ac) + (numel(ap)-1)*var(ap)) / (numel(ac)+numel(ap)-2)); dd = (mean(ac) - mean(ap)) / pooled;
                    cd = (sum(ac(:) > ap(:)') - sum(ac(:) < ap(:)')) / (numel(ac)*numel(ap));
                    rows(end+1,:) = {meths{m}, conds{c}, regs{r}, numel(ac), numel(ap), mean(ac), mean(ap), dd, cd, pv}; %#ok<AGROW>
                end
            end
        end
        GD = cell2table(rows, 'VariableNames', {'Method','Condition','Region','nC','nP','MeanAmp_C','MeanAmp_P','CohenD','CliffDelta','RanksumP'});
        writetable(GD, fullfile(outdir, 'p300_group_difference.csv'));
        fprintf('p300: %d participants (C=%d, P=%d), %.0f s\n', numel(ok), numel(Cs), numel(Ps), toc(t0));

    otherwise
        error('dataset must be dementia, stroke, stroke_saved, AutomaticArtefactRemoval or p300');
end
