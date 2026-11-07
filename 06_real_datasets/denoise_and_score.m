function R = denoise_and_score(X, fs, pairs, id, mode, outs)
% Apply denoising methods to one recording (channels x time, uV) with the parameters of the packaged pipeline
% (WT coif4 hard level by sampling rate; ASR window overlap 0, step 32; EMD 2 IMF) and compute the preservation metrics.
%   mode 'sens' : default settings + sensitivity variants (ASR k, EMD adaptive/3 IMF, WT level and rule)
%   mode 'outs' : evaluate the denoised signals given in struct outs (field = method name)
cfg_path('setup');
X = double(X);
if nargin < 6, outs = struct(); end
out = struct(); out.None = X;
if strcmp(mode, 'outs')
    fn = fieldnames(outs); for q = 1:numel(fn), out.(fn{q}) = double(outs.(fn{q})); end
else
    out.WT = wavelet_denoising_default(X, fs, false, true);
    for L = [5 9], out.(sprintf('WT_L%d', L)) = wt_lvl(X, L, 'Hard'); end
    out.WT_soft = wt_lvl(X, 7, 'Soft');
    nch = size(X,1); E2 = zeros(size(X)); E3 = E2; Ead = E2;
    for ch = 1:nch
        imf = emd(X(ch,:), 'Display', 0); if size(imf,2) < 3, imf(:, end+1:3) = 0; end
        E2(ch,:) = (imf(:,1) + imf(:,2))'; E3(ch,:) = sum(imf(:,1:3), 2)';
        mf = meanfreq(imf, fs); keep = mf >= 1 & mf <= 40; if ~any(keep), keep(1) = true; end
        Ead(ch,:) = sum(imf(:,keep), 2)';
    end
    out.EMD = E2; out.EMD_3IMF = E3; out.EMD_adapt = Ead; out.EMD_WT = wavelet_denoising_default(E2, fs, false, true);
    for k = [5 10 20 30]
        nm = sprintf('ASR_k%d', k);
        try, Y = asr_apply(X, fs, k); catch, Y = nan(size(X)); end
        out.(nm) = Y;
        if k == 5, if all(isfinite(Y(:))), out.ASR_WT = wavelet_denoising_default(Y, fs, false, true); else, out.ASR_WT = Y; end, end
    end
end
pre = metrics_preservation('pre', X, fs, pairs); names = fieldnames(out); rows = cell(numel(names), 1);
for i = 1:numel(names), rows{i} = metrics_preservation('row', id, names{i}, X, out.(names{i}), fs, pairs, pre); end
R = vertcat(rows{:});
end

function Y = wt_lvl(X, lvl, rule)
Y = zeros(size(X));
for ch = 1:size(X,1)
    Y(ch,:) = X(ch,:) - wdenoise(X(ch,:), lvl, 'Wavelet', 'coif4', 'DenoisingMethod', 'Bayes', 'ThresholdRule', rule, 'NoiseEstimate', 'LevelDependent');
end
end
