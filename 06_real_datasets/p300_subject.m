function [R, ERP] = p300_subject(csvfile, id)
% P300 after denoising, Dementia (Muse, 4 channels, 256 Hz).
% markers in the last CSV column (1 Neutral, 2 Positive, 3 Negative), epochs -0.2..0.8 s, baseline -0.2..0, reject > 100 uV.
% Methods are applied to the continuous recording, then the epochs are cut.
cfg_path('setup');
fs = 256;
A = readmatrix(csvfile, 'NumHeaderLines', 1);
X = A(:, 2:5)'; mkr = A(:, 16)'; X(isnan(X)) = 0; mkr(isnan(mkr)) = 0;
Xf = common_filter(X, fs);                              % notch 50 Hz + high-pass 1 Hz, as in the other analyses
% --- methods on the continuous signal
out = struct(); out.None = Xf;
out.WT_repo = wavelet_denoising_default(Xf, fs, false, true);
E2 = zeros(size(Xf)); Ead = E2;
for ch = 1:4
    imf = emd(Xf(ch,:), 'Display', 0); if size(imf,2) < 2, imf(:,2) = 0; end
    E2(ch,:) = (imf(:,1) + imf(:,2))';
    mf = meanfreq(imf, fs); keep = mf >= 1 & mf <= 40; if ~any(keep), keep(1) = true; end
    Ead(ch,:) = sum(imf(:,keep),2)';
end
out.EMD_2IMF = E2; out.EMD_adapt = Ead; out.EMD2_WT = wavelet_denoising_default(E2, fs, false, true);
Y5 = asr_apply(Xf, fs, 5); out.ASR_k5 = Y5; out.ASR5_WT = wavelet_denoising_default(Y5, fs, false, true);
out.ASR_k20 = asr_apply(Xf, fs, 20);
% --- 30 Hz low-pass for ERP analysis, same for every method
[bl, al] = butter(4, 30/(fs/2), 'low'); names = fieldnames(out);
for i = 1:numel(names), out.(names{i}) = filtfilt(bl, al, out.(names{i})')'; end
% --- events: rising edge of the marker channel
on = find(mkr(2:end) ~= 0 & mkr(1:end-1) == 0) + 1; code = mkr(on);
pre = round(0.2*fs); post = round(0.8*fs); t = (-pre:post)/fs*1000;
ok = on - pre >= 1 & on + post <= size(Xf,2) & ismember(code, [1 2 3]); on = on(ok); code = code(ok);
Ref = out.None; keep = false(numel(on),1);
for e = 1:numel(on)
    ep = Ref(:, on(e)-pre:on(e)+post); ep = ep - mean(ep(:, 1:pre), 2);
    keep(e) = all(abs(ep(:)) <= 100);                  % same trial set for every method (criterion on the unprocessed epoch)
end
on = on(keep); code = code(keep);
conds = {'Neutral', 'Positive', 'Negative'}; regions = {'Frontal', 'Temporal'}; rch = {[2 3], [1 4]};   % AF7/AF8 ; TP9/TP10
w = t >= 300 & t <= 550; bwin = t < 0;
rows = {}; ERP = struct('t', t);
for i = 1:numel(names)
    Y = out.(names{i});
    for c = 1:3
        idx = on(code == c); if numel(idx) < 5, continue; end
        for r = 1:2
            E = zeros(numel(idx), numel(t));
            for e = 1:numel(idx)
                seg = mean(Y(rch{r}, idx(e)-pre:idx(e)+post), 1); E(e,:) = seg - mean(seg(bwin));
            end
            avg = mean(E, 1); ERP.(names{i}).(conds{c}).(regions{r}) = avg; [pk, ip] = max(avg(w)); tw = t(w); lat = tw(ip);
            snr = pk / (std(avg(bwin)) + eps); var_tr = mean(std(E(:, w), 0, 1)); ecorr = corr(avg(:), ERP.None.(conds{c}).(regions{r})(:));
            rows(end+1,:) = {id, names{i}, conds{c}, regions{r}, numel(idx), pk, lat, snr, var_tr, mean(avg(w)), ecorr}; %#ok<AGROW>
        end
    end
end
R = cell2table(rows, 'VariableNames', {'Rec','Method','Condition','Region','nTrials','P3_peak_uV','P3_latency_ms','ERP_SNR','TrialSD_uV','P3_meanAmp_uV','ERP_corr'});
end



