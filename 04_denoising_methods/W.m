%% W: wavelet denoising (WT), as used in the paper
%
% Input : EEG_original_int.data  (channels x samples; after common_filter, and after interpolation of bad
%                                 channels when the dataset script does it), EEG.srate, EEG.nbchan
% Output: EEG.data_W, EEG.residual_W
%
% Discrete wavelet transform (MATLAB wdenoise) per channel: Coiflet 4, decomposition level 8 / 7 / 6 for sampling rates
% above 500 Hz / between 250 and 500 Hz / up to 250 Hz, level-dependent Bayesian threshold with the hard rule. The
% component removed by the threshold (the artefact) is subtracted from the signal (see 04_denoising_methods/wavelet_denoising_default.m).

EEG_int.data_filt = EEG_original_int.data;
denoised_EEG = zeros(size(EEG_int.data_filt));
residual_noise = zeros(size(EEG_int.data_filt));

srate = EEG.srate;
paradigmERP = false;
softThresh = true;      % the soft rule is applied only if paradigmERP is true; otherwise the hard rule is used

for ch = 1:EEG.nbchan
    denoised_EEG(ch, :) = wavelet_denoising_default(EEG_int.data_filt(ch, :), srate, paradigmERP, softThresh);
    residual_noise(ch, :) = EEG_int.data_filt(ch, :) - denoised_EEG(ch, :);
end

EEG.data_W = denoised_EEG;
EEG.residual_W = residual_noise;
