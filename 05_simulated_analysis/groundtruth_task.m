function [T, info, timing] = groundtruth_task(seed, scale, layer)
% One simulated subject at one artefact amplitude. layer: 'full' (all variants) or 'main' (default settings only).
% Metrics are vs the CLEAN ground truth. *_aff = averaged over contaminated channels only (comparable across densities).
cfg_path('setup');
fs = 256;
[clean, contam, info] = sim_make(seed, scale);
Cf = common_filter(clean, fs); Xf = common_filter(contam, fs);
Channels = info.channels;
ix = @(names) find(ismember(Channels, names));
ex16 = {'AFz','Fz','FT7','FC3','FCz','FC4','FT8','Cz','TP7','CP3','CPz','CP4','TP8','Pz','OZ'};
c16 = Channels(~ismember(Channels, ex16));
c8  = c16(~ismember(c16, {'F7','F3','F4','F8','P7','P3','P4','P8'}));
c4  = c8(~ismember(c8, {'C3','C4','O1','O2'}));
cfg = struct('n', {31, 16, 8, 4}, 'rows', {1:31, ix(c16), ix(c8), ix(c4)});
full = strcmp(layer, 'full');
timing = struct();

out = struct();
out.None = Xf;
tic; out.WT_repo = wavelet_denoising_default(Xf, fs, false, true); timing.WT_repo = toc;
if full
    tic;
    out.WT_repo_soft = wt_apply(Xf, 'bior4.4', 7, 'Soft', 'resid');
    out.WT_repo_L5   = wt_apply(Xf, 'bior4.4', 5, 'Hard', 'resid');
    out.WT_repo_L9   = wt_apply(Xf, 'bior4.4', 9, 'Hard', 'resid');
    out.WT_repo_db4  = wt_apply(Xf, 'db4',   7, 'Hard', 'resid');
    out.WT_repo_bior44 = wt_apply(Xf, 'bior4.4', 7, 'Hard', 'resid');
    out.WT_repo_sym4 = wt_apply(Xf, 'sym4',  7, 'Hard', 'resid');
    out.WT_repo_univ = wt_apply(Xf, 'bior4.4', 7, 'Hard', 'resid', 'UniversalThreshold');
    timing.WT_variants = toc;
end
tic;
out.WT_den_hard = wt_apply(Xf, 'bior4.4', 7, 'Hard', 'den');
out.WT_den_soft = wt_apply(Xf, 'bior4.4', 7, 'Soft', 'den');
timing.WT_den = toc;
tic;
nch = size(Xf,1); E2 = zeros(size(Xf)); E3 = E2; E4 = E2; Erm2 = E2; Ead = E2; nkeep = zeros(nch,1);
for ch = 1:nch
    x = Xf(ch,:); imf = emd(x, 'Display', 0);
    if size(imf,2) < 4, imf(:, end+1:4) = 0; end
    E2(ch,:) = (imf(:,1) + imf(:,2))'; E3(ch,:) = sum(imf(:,1:3),2)'; E4(ch,:) = sum(imf(:,1:4),2)';
    Erm2(ch,:) = x - E2(ch,:);
    mf = meanfreq(imf, fs); keep = mf >= 1 & mf <= 40; if ~any(keep), keep(1) = true; end
    Ead(ch,:) = sum(imf(:,keep),2)'; nkeep(ch) = sum(keep);
end
out.EMD_2IMF = E2;
if full, out.EMD_3IMF = E3; out.EMD_4IMF = E4; out.EMD_rm2 = Erm2; end
out.EMD_adapt = Ead; info.emd_adapt_kept_mean = mean(nkeep);
timing.EMD = toc;
tic; out.EMD2_WT = wavelet_denoising_default(E2, fs, false, true); timing.EMD2_WT = toc;
if full
    Wout = out.WT_repo; WE = zeros(size(Wout));
    for ch = 1:nch, imf = emd(Wout(ch,:), 'Display', 0); WE(ch,:) = (imf(:,1) + imf(:,2))'; end
    out.WT_EMD2 = WE;
end

rowsT = {};
names = fieldnames(out);
for m = 1:numel(names)
    for c = 1:numel(cfg)
        rowsT(end+1,:) = mk(seed, scale, cfg(c), names{m}, Cf, out.(names{m}), Xf, fs, info.affected); %#ok<AGROW>
    end
end
if full, ks = [5 10 20 30]; else, ks = [5 20]; end
for c = 1:numel(cfg)
    Xs = Xf(cfg(c).rows,:);
    for k = ks
        tic; Y = asr_apply(Xs, fs, k); timing.(sprintf('ASR_k%d_ch%d', k, cfg(c).n)) = toc;
        Yfull = Xf; Yfull(cfg(c).rows,:) = Y;
        rowsT(end+1,:) = mk(seed, scale, cfg(c), sprintf('ASR_k%d', k), Cf, Yfull, Xf, fs, info.affected); %#ok<AGROW>
        if k == 5
            Yfull(cfg(c).rows,:) = wavelet_denoising_default(Y, fs, false, true);
            rowsT(end+1,:) = mk(seed, scale, cfg(c), 'ASR5_WT', Cf, Yfull, Xf, fs, info.affected); %#ok<AGROW>
        end
    end
end
T = cell2table(rowsT, 'VariableNames', {'Subject','Scale','Nch','Method','SNR','RRMSE_t','RRMSE_f','CC', ...
    'SNR_aff','RRMSE_t_aff','RRMSE_f_aff','CC_aff','SNR_repo','RRMSE_t_repo','CC_repo'});
end

function r = mk(seed, scale, cfg, name, Cf, D, Xf, fs, affected)
M = metrics_groundtruth(Cf, D, Xf, fs, cfg.rows);
a = ismember(cfg.rows, affected);
r = {seed, scale, cfg.n, name, mean(M.SNR), mean(M.RRMSE_t), mean(M.RRMSE_f), mean(M.CC), ...
     mean(M.SNR(a)), mean(M.RRMSE_t(a)), mean(M.RRMSE_f(a)), mean(M.CC(a)), ...
     mean(M.SNR_r(a)), mean(M.RRMSE_t_r(a)), mean(M.CC_r(a))};
end


