function varargout = metrics_preservation(action, varargin)
% Reference-free preservation metrics for real recordings.
%   pre = metrics_preservation('pre', X, fs, pairs)                    % statistics of the unprocessed signal
%   R   = metrics_preservation('row', id, name, X, Y, fs, pairs, pre)  % one table row for method "name"
% Metrics (post vs pre): GA at +-50/100/150 uV, band power ratio (delta, theta, alpha, beta), RMS ratio (sqrt of variance ratio),
% kurtosis, spectral entropy, preservation of the inter-channel correlation structure, preservation of the autocorrelation
% function (lags 0-1 s), and pdBSI per band (pair-derived Brain Symmetry Index,
% frequency by frequency |R-L|/(R+L| averaged within the band and over homologous pairs).
switch action
    case 'pre', varargout{1} = stats(varargin{1}, varargin{2}, varargin{3});
    case 'row', varargout{1} = row(varargin{:});
end
end

function B = stats(Y, fs, pairs)
nfft = 2*round(fs); [P, f] = pwelch(Y', hann(nfft), nfft/2, nfft, fs);   % freq x channels
bands = [1 4; 4 8; 8 13; 13 30];
B.bp = zeros(size(Y,1), 4);
for b = 1:4, m = f >= bands(b,1) & f < bands(b,2); B.bp(:,b) = trapz(f(m), P(m,:), 1)'; end
m = f >= 1 & f <= 45; Pm = P(m,:); p = Pm ./ sum(Pm, 1); N = size(Pm,1);
B.entropy = (-sum(p .* log2(p + eps), 1) / log2(N))';
B.pdbsi = nan(1, 4);
for b = 1:4
    mm = f >= bands(b,1) & f < bands(b,2); v = nan(size(pairs,1), 1);
    for j = 1:size(pairs,1)
        l = P(mm, pairs(j,1)); r = P(mm, pairs(j,2)); v(j) = mean(abs(r - l) ./ (r + l + eps));
    end
    B.pdbsi(b) = mean(v, 'omitnan');
end
L = round(fs); B.acf = zeros(size(Y,1), L);
for ch = 1:size(Y,1)
    a = xcorr(Y(ch,:) - mean(Y(ch,:)), L, 'coeff'); B.acf(ch,:) = a(L+2:end);
end
end

function ga = ga_score(Y, fs, T)
W = round(2*fs); K = floor(size(Y,2)/W); edges = -1000:10:1000; sc = zeros(K,1);
nb0 = find(edges >= -T, 1, 'first'); nb1 = find(edges <= T, 1, 'last') - 1;
for k = 1:K
    w = Y(:, (k-1)*W+1:k*W); c = histcounts(w(:), edges); tot = sum(c);
    if tot == 0, sc(k) = 0; continue; end
    s1 = sum(c(nb0:nb1)) / tot; s2 = sum(c(nb0:nb1) > 0) / sum(c > 0); sc(k) = (s1 + s2) / 2 * 100;
end
ga = mean(sc);
end

function r = row(id, name, X, Y, fs, pairs, pre)
nanrow = nan(1, 21);
if any(~isfinite(Y(:))) || size(Y,2) ~= size(X,2), vals = nanrow; else
    post = stats(Y, fs, pairs);
    bpr = mean(post.bp ./ (pre.bp + eps), 1);
    rmsr = mean(sqrt(mean(Y.^2,2)) ./ (sqrt(mean(X.^2,2)) + eps));
    kurt = mean(kurtosis(Y, 1, 2));
    C1 = corrcoef(X'); C2 = corrcoef(Y'); iu = find(triu(true(size(C1)), 1)); a = C1(iu); b = C2(iu); ok = isfinite(a) & isfinite(b);
    if sum(ok) >= 2 && std(b(ok)) > 0, covr = corr(a(ok), b(ok)); else, covr = 0; end
    ac = zeros(size(Y,1),1);
    for ch = 1:size(Y,1), u = pre.acf(ch,:)'; v = post.acf(ch,:)'; if std(v) > 0 && std(u) > 0, ac(ch) = corr(u, v); end, end
    vals = [ga_score(Y,fs,50) ga_score(Y,fs,100) ga_score(Y,fs,150) bpr rmsr kurt mean(post.entropy) mean(pre.entropy) covr mean(ac) ...
            post.pdbsi pre.pdbsi];
    
end
r = table(string(id), string(name), 'VariableNames', {'Rec','Method'});
vn = {'GA50','GA100','GA150','BPR_delta','BPR_theta','BPR_alpha','BPR_beta','RMS_ratio','Kurtosis','Entropy','Entropy_pre','CovPreserve','ACF_corr', ...
      'pdBSI_delta','pdBSI_theta','pdBSI_alpha','pdBSI_beta','pdBSI_pre_delta','pdBSI_pre_theta','pdBSI_pre_alpha','pdBSI_pre_beta'};
for i = 1:numel(vn), r.(vn{i}) = vals(i); end
end

