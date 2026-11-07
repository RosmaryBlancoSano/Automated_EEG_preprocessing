function R = score_auto_subject(setfile, id)
% One AutomaticArtefactRemoval participant (whole recording). Channel-wise methods are computed once on the 27 EEG channels;
% ASR (which depends on the channel subset) is computed per configuration.
cfg_path('setup');
[X, fs, Einfo] = read_eeglab_set(setfile); labels = Einfo.labels;
keep = ~ismember(labels, {'EOG1','EOG2','EOG3'}); X = X(keep,:); labels = labels(keep);
X = common_filter(X, fs);                                         % notch 50 Hz + high-pass 1 Hz
% channel-wise outputs (27 channels)
cw = struct();
cw.WT = wavelet_denoising_default(X, fs, false, true);
cw.WT_L5 = wt_lvl(X, 5, 'Hard'); cw.WT_L9 = wt_lvl(X, 9, 'Hard'); cw.WT_soft = wt_lvl(X, 7, 'Soft');
E2 = zeros(size(X)); E3 = E2; Ead = E2;
for ch = 1:size(X,1)
    imf = emd(X(ch,:), 'Display', 0); if size(imf,2) < 3, imf(:, end+1:3) = 0; end
    E2(ch,:) = (imf(:,1) + imf(:,2))'; E3(ch,:) = sum(imf(:,1:3), 2)';
    mf = meanfreq(imf, fs); k = mf >= 1 & mf <= 40; if ~any(k), k(1) = true; end
    Ead(ch,:) = sum(imf(:,k), 2)';
end
cw.EMD = E2; cw.EMD_3IMF = E3; cw.EMD_adapt = Ead; cw.EMD_WT = wavelet_denoising_default(E2, fs, false, true);
cfgs = {27, labels; 16, {'Fp1','Fp2','F7','F3','F4','F8','T7','C3','C4','T8','P7','P3','P4','P8','O1','O2'}; 8, {'Fp1','Fp2','T7','C3','C4','T8','O1','O2'}; 4, {'F3','F4','T7','T8'}};
allpairs = {'Fp1','Fp2'; 'F7','F8'; 'F3','F4'; 'T7','T8'; 'C3','C4'; 'P7','P8'; 'P3','P4'; 'O1','O2'; 'FC5','FC6'; 'FC1','FC2'; 'CP5','CP6'; 'CP1','CP2'};
rows = {};
for c = 1:size(cfgs,1)
    idx = find(ismember(labels, cfgs{c,2})); names = labels(idx); Xs = X(idx,:);
    pairs = [];
    for q = 1:size(allpairs,1), a = find(strcmp(names, allpairs{q,1})); b = find(strcmp(names, allpairs{q,2})); if ~isempty(a) && ~isempty(b), pairs(end+1,:) = [a b]; end, end %#ok<AGROW>
    outs = struct(); fn = fieldnames(cw); for q = 1:numel(fn), outs.(fn{q}) = cw.(fn{q})(idx,:); end
    for k = [5 10 20 30]
        try, Y = asr_apply(Xs, fs, k); catch, Y = nan(size(Xs)); end
        outs.(sprintf('ASR_k%d', k)) = Y;
        if k == 5, if all(isfinite(Y(:))), outs.ASR_WT = wavelet_denoising_default(Y, fs, false, true); else, outs.ASR_WT = Y; end, end
    end
    R1 = denoise_and_score(Xs, fs, pairs, id, 'outs', outs); R1.Nch = repmat(cfgs{c,1}, height(R1), 1); rows{end+1} = R1; %#ok<AGROW>
end
R = vertcat(rows{:});
end

function Y = wt_lvl(X, lvl, rule)
Y = zeros(size(X));
for ch = 1:size(X,1)
    Y(ch,:) = X(ch,:) - wdenoise(X(ch,:), lvl, 'Wavelet', 'coif4', 'DenoisingMethod', 'Bayes', 'ThresholdRule', rule, 'NoiseEstimate', 'LevelDependent');
end
end
