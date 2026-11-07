function Y = asr_apply(X, fs, k)
%ASR_APPLY  Artefact subspace reconstruction as used in the paper (single implementation, called by ASR.m and by the evaluation scripts).
%   Y = asr_apply(X, fs, k)   X: channels x samples (after common_filter), k: cutoff (the paper uses 5 as default).
%   ASR is calibrated on the recording itself. Parameters: blocksize 10, windows of 0.5 s without overlap, maximum dropout fraction 0.1,
%   minimum clean fraction 0.25, statistics updated every 32 samples, lookahead 0.25 s, maximum dimensionality 0.66. The maximum memory
%   is passed explicitly (256 MB) so that the Java memory query of clean_rawdata is not needed. The output is delayed by the processing
%   window and is realigned with the input by cross-correlation (lag searched within +-2 s).
cfg_path('setup');   % asr_calibrate and asr_process are in 04_denoising_methods/
blocksize = 10; window_len = 0.5; window_overlap = 0.0; maxdrop = 0.1; mincl = 0.25;
state = asr_calibrate(X, fs, k, blocksize, [], [], window_len, window_overlap, maxdrop, mincl, 256);
lookahead = window_len/2; stepsize = 32; maxdims = 0.66;
evalc('[Y, ~] = asr_process(X, fs, state, window_len, lookahead, stepsize, maxdims, 256);');
% realign with the input (cross-correlation on channel 1)
[xc, lags] = xcorr(X(1,:), Y(1,:), 2*fs, 'coeff'); [~, im] = max(xc);
Y = circshift(Y, lags(im), 2);
end