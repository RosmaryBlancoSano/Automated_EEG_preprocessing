function F = stroke_features(X, fs)
% Stroke biomarkers, computed on non-overlapping 10-s windows and averaged over windows:
% relative band power (1-4, 4-8, 8-13, 13-30 Hz; relative to 1-30 Hz), Delta/Alpha, Theta/Alpha, BSI per band,
% RMS, kurtosis. PSD: Hamming-windowed periodogram over the whole window (n_per_seg = window), zero-padded FFT.
% Channels: FC3 FC4 CP3 CP4 FT7 FT8 TP7 TP8; homologous pairs (1,2) (3,4) (5,6) (7,8).
bands = [1 4; 4 8; 8 13; 13 30]; pairs = [1 2; 3 4; 5 6; 7 8];
W = round(10*fs); K = floor(size(X,2)/W);
F = struct('RelDelta', NaN, 'RelTheta', NaN, 'RelAlpha', NaN, 'RelBeta', NaN, 'DeltaAlpha', NaN, 'ThetaAlpha', NaN, ...
           'BSI_Delta', NaN, 'BSI_Theta', NaN, 'BSI_Alpha', NaN, 'BSI_Beta', NaN, 'RMS', NaN, 'Kurt', NaN, 'nWin', K);
if K < 1, return; end
nfft = 2^nextpow2(W); hw = 0.54 - 0.46*cos(2*pi*(0:W-1)'/(W-1)); f = (0:nfft/2)'*fs/nfft;
acc = nan(K, 12);
for k = 1:K
    seg = X(:, (k-1)*W+1:k*W)';                                  % W x ch
    P = abs(fft((seg - mean(seg,1)) .* hw, nfft)).^2; P = P(1:nfft/2+1, :);
    bp = zeros(4, size(seg,2));
    for b = 1:4, m = f >= bands(b,1) & f < bands(b,2); bp(b,:) = trapz(f(m), P(m,:), 1); end
    tot = sum(bp, 1); rel = bp ./ (tot + eps);
    bsi = zeros(1,4);
    for b = 1:4
        m = f >= bands(b,1) & f < bands(b,2); s = sum(P(m,:), 1); v = [];
        for j = 1:size(pairs,1), l = s(pairs(j,1)); r = s(pairs(j,2)); if l + r > 0, v(end+1) = abs((l-r)/(l+r)); end, end %#ok<AGROW>
        bsi(b) = mean(v);
    end
    acc(k,:) = [mean(rel(1,:)) mean(rel(2,:)) mean(rel(3,:)) mean(rel(4,:)) mean(rel(1,:) ./ (rel(3,:)+eps)) mean(rel(2,:) ./ (rel(3,:)+eps)) bsi ...
                mean(sqrt(mean(seg.^2,1))) mean(kurtosis(seg,1,1))];
end
v = mean(acc, 1, 'omitnan');
F.RelDelta = v(1); F.RelTheta = v(2); F.RelAlpha = v(3); F.RelBeta = v(4); F.DeltaAlpha = v(5); F.ThetaAlpha = v(6);
F.BSI_Delta = v(7); F.BSI_Theta = v(8); F.BSI_Alpha = v(9); F.BSI_Beta = v(10); F.RMS = v(11); F.Kurt = v(12);
end
