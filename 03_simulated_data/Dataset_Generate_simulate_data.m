%% ---------------------------------------------------------------------
%% Dataset: semi-simulated EEG with ground truth
%% ---------------------------------------------------------------------
% For each subject (random seed = subject number, so the data can be regenerated exactly):
%   1. clean resting-state EEG, 31 channels, generated with fake_EEG.m (after Yeung et al.): a 10 Hz alpha rhythm over posterior
%      channels and a 20 Hz beta rhythm over fronto-central channels, jittered over 973 trials of 100 samples, projected with a
%      dipole-based spatial distribution, plus noise with the power spectrum of human EEG; 97300 samples, stored at 256 Hz (380 s)
%   2. 100 EOG and 100 EMG epochs of EEGdenoiseNet (https://github.com/ncclabsustech/EEGdenoiseNet), concatenated, resampled and
%      added linearly: EOG to frontal channels, EMG to temporal and central channels (sim_make.m)
%   3. channel subsets of 16, 8 and 4 channels (nested subsets of the 10-20 montage)
% The CLEAN signal is saved too (_clean files): it is the ground truth of Simulated_groundtruth_evaluation.m.
% Output: results/SimulatedEEG/subjectN/EEG_<n>_subjectN.mat (contaminated) and EEG_<n>_subjectN_clean.mat

clear; clc;
addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');
base_path = fullfile(cfg_path('results'), 'SimulatedEEG');
num_subjects = 10;
amplitude = 1;        % multiplier of the artefact amplitude (1 = original)
fs = 256;

exclude_16 = {'AFz', 'Fz', 'FT7', 'FC3', 'FCz', 'FC4', 'FT8', 'Cz', 'TP7', 'CP3', 'CPz', 'CP4', 'TP8', 'Pz', 'OZ'};
exclude_8  = {'F7', 'F3', 'F4', 'F8', 'P7', 'P3', 'P4', 'P8'};
exclude_4  = {'C3', 'C4', 'O1', 'O2'};

for subj = 1:num_subjects
    [clean, contam, info] = sim_make(subj, amplitude);
    Channels = info.channels;
    Channels_16 = Channels(~ismember(Channels, exclude_16));
    Channels_8  = Channels_16(~ismember(Channels_16, exclude_8));
    Channels_4  = Channels_8(~ismember(Channels_8, exclude_4));
    sets = {31, Channels; 16, Channels_16; 8, Channels_8; 4, Channels_4};

    subj_folder = fullfile(base_path, sprintf('subject%d', subj));
    if ~exist(subj_folder, 'dir'), mkdir(subj_folder); end
    for c = 1:size(sets, 1)
        keep = ismember(Channels, sets{c, 2});
        for kind = {'', '_clean'}
            X = contam; if strcmp(kind{1}, '_clean'), X = clean; end
            EEGs.data = X(keep, :);
            EEGs.chanlocs = struct('labels', sets{c, 2});
            EEGs.srate = fs;
            EEGs.trials = 1;
            EEGs.event = 0;
            EEGs.nbchan = sets{c, 1};
            save(fullfile(subj_folder, sprintf('EEG_%d_subject%d%s.mat', sets{c, 1}, subj, kind{1})), 'EEGs');
        end
    end
    fprintf('Subject %d: signal-to-artefact ratio of the contaminated channels %.1f dB (median)\n', subj, median(info.SAR_dB));
end
