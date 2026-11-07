%% ASR: artefact subspace reconstruction (clean_rawdata), as used in the paper
%
% Input : EEG_original_int.data  (channels x samples, after common_filter), EEG.srate
%         asr_cutoff (optional)  cutoff k; default 5 (the most aggressive value usually suggested)
% Output: EEG.data_ASR, EEG.residual_ASR
%
% ASR is calibrated on the recording itself (no separate clean calibration data). The parameters (windows of 0.5 s without overlap,
% statistics updated every 32 samples, lookahead 0.25 s, maximum dropout fraction 0.1, minimum clean fraction 0.25, maximum
% dimensionality 0.66) and the realignment of the output with the input by cross-correlation are in 04_denoising_methods/asr_apply.m, which is
% the only implementation of ASR in this repository (also used by Simulated_groundtruth_evaluation and Real_datasets_evaluation).

EEG_int.data_filt = EEG_original_int.data;
if ~exist('asr_cutoff', 'var') || isempty(asr_cutoff), asr_cutoff = 5; end

denoised_EEG = asr_apply(EEG_int.data_filt, EEG.srate, asr_cutoff);

EEG.data_ASR = denoised_EEG;
EEG.residual_ASR = EEG_int.data_filt - denoised_EEG;