%% EMD: empirical mode decomposition, "2IMF" rule, as used in the paper
%
% Input : EEG_original_int.data  (channels x samples, after common_filter), EEG.srate, EEG.nbchan
% Output: EEG.data_EMD, EEG.residual_EMD
%
% Each channel is decomposed into intrinsic mode functions (IMFs) and reconstructed from the first two, which are the
% fastest oscillations (04_denoising_methods/EMD_denoise_EEG_core.m, after Mert and Akan). The slower IMFs, which contain the ocular
% artefact, are discarded: slow neural activity is removed with them.

EEG_int.data_filt = EEG_original_int.data;
fs = EEG.srate;
denoised_EEG = zeros(size(EEG_int.data_filt));
residual_noise = zeros(size(EEG_int.data_filt));

noise_type = 'eog';         % type of noise ('emg' or 'eog')
classify_method = '2IMF';   % noise classification method ('threshold' or '2IMF')

for ch = 1:EEG.nbchan
    denoised_EEG(ch, :) = EMD_denoise_EEG_core(EEG_int.data_filt(ch, :), fs, noise_type, classify_method);
    residual_noise(ch, :) = EEG_int.data_filt(ch, :) - denoised_EEG(ch, :);
end

EEG.data_EMD = denoised_EEG;
EEG.residual_EMD = residual_noise;
