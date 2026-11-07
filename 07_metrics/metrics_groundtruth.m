function M = metrics_groundtruth(C, D, X, fs, rows)
% Per-channel metrics, averaged over 2 s windows with 50% overlap.
%   C: clean ground truth, D: denoised output, X: contaminated input (all channels x time)
%   Ground-truth metrics (vs C):   SNR (dB), RRMSE_t, RRMSE_f, CC
%   "Repo-style" metrics (vs X):   SNR_r, RRMSE_t_r, CC_r (what compareDenoisingMethods.m computes)
% rows: channel indices to evaluate
W = round(2*fs); step = W/2; T = size(C,2);
K = floor((T-W)/step) + 1;
idx = (1:W)' + (0:K-1)*step;
hw = 0.5 - 0.5*cos(2*pi*(0:W-1)'/(W-1));
nb = floor(W/2) + 1;
nr = numel(rows);
M.SNR = zeros(1,nr); M.RRMSE_t = M.SNR; M.RRMSE_f = M.SNR; M.CC = M.SNR;
M.SNR_r = M.SNR; M.RRMSE_t_r = M.SNR; M.CC_r = M.SNR;
for i = 1:nr
    ch = rows(i);
    c = reshape(C(ch, idx(:)), W, K); d = reshape(D(ch, idx(:)), W, K); x = reshape(X(ch, idx(:)), W, K);
    [M.SNR(i), M.RRMSE_t(i), M.RRMSE_f(i), M.CC(i)] = win_metrics(c, d, hw, nb, true);
    [M.SNR_r(i), M.RRMSE_t_r(i), ~, M.CC_r(i)] = win_metrics(x, d, hw, nb, false);
end
end

function [snr, rrt, rrf, cc] = win_metrics(c, d, hw, nb, power_snr)
e = d - c;
rc = sqrt(mean(c.^2, 1)); re = sqrt(mean(e.^2, 1));
if power_snr
    snrw = repmat(10*log10(mean(rc.^2) / (mean(re.^2) + eps)), size(rc));   % global output SNR vs ground truth (no per-window clipping artefacts)
else
    snrw = 10*log10(rc ./ (re + eps));                    % as in compareDenoisingMethods.m
end
rrtw = re ./ (rc + eps);
Pc = abs(fft(c .* hw)).^2; Pd = abs(fft(d .* hw)).^2; Pc = Pc(1:nb,:); Pd = Pd(1:nb,:);
rrfw = sqrt(mean((Pd - Pc).^2, 1)) ./ (sqrt(mean(Pc.^2, 1)) + eps);
zc = (c - mean(c,1)); zd = (d - mean(d,1));
den = sqrt(sum(zc.^2,1) .* sum(zd.^2,1));
ccw = sum(zc .* zd, 1) ./ den; ccw(den == 0 | ~isfinite(ccw)) = 0;
snrw = min(max(snrw,-60),60); snr = mean(snrw); rrt = mean(rrtw); rrf = mean(rrfw); cc = mean(ccw);
end


