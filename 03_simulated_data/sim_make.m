function [clean, contam, info] = sim_make(seed, scale)
% Regenerates a semi-simulated EEG subject keeping the CLEAN ground truth.
% Same recipe as Dataset_Generate_simulate_data.m (fake_EEG + EOG/EMG from EEGdenoiseNet,
% linear mixing), with an extra artefact amplitude multiplier "scale" (1 = original).
cfg_path('setup'); base = cfg_path('eegdenoisenet');

persistent EMG EOG
if isempty(EMG)
    get_eegdenoisenet(base);   % downloads the public EEGdenoiseNet epochs the first time
    EMG = load(fullfile(base, 'EMG_all_epochs.mat')); EMG = EMG.EMG_all_epochs;
    EOG = load(fullfile(base, 'EOG_all_epochs.mat')); EOG = EOG.EOG_all_epochs;
end
rng(seed, 'twister');
clean = fake_EEG;            % 31 x 97300, silence the per-channel printouts
Channels = {'Fp1','Fp2','AFz','F7','F3','Fz','F4','F8','FT7','FC3','FCz','FC4','FT8','T7','C3','Cz','C4','T8','TP7','CP3','CPz','CP4','TP8','P7','P3','Pz','P4','P8','O1','OZ','O2'};
A = randsample(size(EMG,1), 100); B = randsample(size(EOG,1), 100);
emg = reshape(EMG(A,:)', 1, []); eog = reshape(EOG(B,:)', 1, []);
n = size(clean, 2);
emg_r = resample(emg', n, numel(emg))'; eog_r = resample(eog', n, numel(eog))';
front_idx = ismember(Channels, {'Fp1','Fp2','AFz','F7','F3','Fz','F4','F8'});
temp_idx  = ismember(Channels, {'TP7','TP8','T7','T8','C3','C4'});
U = [1.2; 1.0; 0.8; 0.9; 1.1; 0.7] / 1500;
P = [1; 1; 0.73; 0.1; 0.3; 0.2; 0.3; 0.1] / 3;
mus = scale * (U * emg_r);            % 6 x n
ocu = scale * (P * eog_r) * -1;       % 8 x n
contam = clean;
contam(front_idx, :) = clean(front_idx, :) + ocu;
contam(temp_idx, :)  = clean(temp_idx, :) + mus;
art = contam - clean;
aff = find(front_idx | temp_idx);
info.channels = Channels;
info.affected = aff;
info.SAR_dB = 10*log10(var(clean(aff,:), 0, 2) ./ var(art(aff,:), 0, 2));
info.art_p99 = prctile(abs(art(aff,:)), 99, 2);
info.p99_clean = prctile(abs(clean(:)), 99);
info.p99_contam = prctile(abs(contam(:)), 99);
info.frac_within100_clean = mean(abs(clean(:)) <= 100);
info.frac_within100_contam = mean(abs(contam(:)) <= 100);
info.scale = scale; info.seed = seed;
end

