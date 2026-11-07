function Y = wt_apply(X, wname, lvl, rule, mode, method)
% Wavelet thresholding per channel with wdenoise.
%   mode 'resid' : output = x - wdenoise(x)   (what wavelet_denoising_default.m does)
%   mode 'den'   : output = wdenoise(x)       (conventional wavelet denoising)
if nargin < 6, method = 'Bayes'; end
Y = zeros(size(X));
for ch = 1:size(X,1)
    x = double(X(ch,:));
    d = wdenoise(x, lvl, 'Wavelet', wname, 'DenoisingMethod', method, 'ThresholdRule', rule, 'NoiseEstimate', 'LevelDependent');
    if strcmp(mode, 'resid'), Y(ch,:) = x - d; else, Y(ch,:) = d; end
end
end
