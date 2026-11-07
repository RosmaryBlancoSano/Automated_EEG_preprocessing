%% ASR_W: ASR followed by wavelet denoising (ASR -> WT), as used in the paper
%
% Input : EEG_original_int.data  (channels x samples, after common_filter), EEG.srate, EEG.nbchan
% Output: EEG.data_ASR_W, EEG.residual_ASR_W
%
% Stage 1: ASR (ASR.m, cutoff 5). Stage 2: the wavelet denoising of W.m applied to the OUTPUT of stage 1.

EEG_int.data_filt = EEG_original_int.data;

ASR                                   % stage 1 -> EEG.data_ASR (realigned with the input)
denoised_EEG_temp = EEG.data_ASR;

srate = EEG.srate;
paradigmERP = false;
softThresh = true;      % the soft rule is applied only if paradigmERP is true; otherwise the hard rule is used

denoised_EEG = zeros(size(EEG_int.data_filt));
residual_noise = zeros(size(EEG_int.data_filt));
for ch = 1:EEG.nbchan
    denoised_EEG(ch, :) = wavelet_denoising_default(denoised_EEG_temp(ch, :), srate, paradigmERP, softThresh);
    residual_noise(ch, :) = EEG_int.data_filt(ch, :) - denoised_EEG(ch, :);
end

EEG.data_ASR_W = denoised_EEG;
EEG.residual_ASR_W = residual_noise;
