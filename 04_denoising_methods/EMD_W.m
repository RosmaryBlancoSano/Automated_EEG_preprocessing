%% EMD_W: EMD ("2IMF") followed by wavelet denoising (EMD -> WT), as used in the paper
%
% Input : EEG_original_int.data  (channels x samples, after common_filter), EEG.srate, EEG.nbchan
% Output: EEG.data_EMD_W, EEG.residual_EMD_W
%
% Stage 1: EMD keeping the first two IMFs. Stage 2: the wavelet denoising of W.m applied to the OUTPUT of stage 1.

EEG_int.data_filt = EEG_original_int.data;
denoised_EEG = zeros(size(EEG_int.data_filt));
residual_noise = zeros(size(EEG_int.data_filt));
denoised_temp = zeros(size(EEG_int.data_filt));

fs = EEG.srate;
paradigmERP = false;
softThresh = true;      % the soft rule is applied only if paradigmERP is true; otherwise the hard rule is used

for ch = 1:EEG.nbchan
    % Stage 1: EMD, first two IMFs
    denoised_temp(ch, :) = EMD_denoise_EEG_core(EEG_int.data_filt(ch, :), fs, 'eog', '2IMF');
    % Stage 2: wavelet denoising of the stage-1 output
    denoised_EEG(ch, :) = wavelet_denoising_default(denoised_temp(ch, :), fs, paradigmERP, softThresh);
    residual_noise(ch, :) = EEG_int.data_filt(ch, :) - denoised_EEG(ch, :);
end

EEG.data_EMD_W = denoised_EEG;
EEG.residual_EMD_W = residual_noise;
