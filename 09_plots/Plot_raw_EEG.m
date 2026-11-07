addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');   % folders of the repository (see cfg_path.m)

clear; clc;

% Set Paths
% EEGLAB (https://sccn.ucsd.edu/eeglab, 2022.1 or later) must be on the MATLAB path and started (eeglab nogui) before running
% EEGLAB (2022.1 or later) must be on the MATLAB path before running
base_path = cfg_path('stroke');

% Subject number
subj_folder = 'sub-0001';  
subj_path = fullfile(base_path, subj_folder);

% Load .set file
eeg_file = dir(fullfile(subj_path, '*Original.set'));

if ~isempty(eeg_file)
    file_name = eeg_file(1).name;
    file_path = fullfile(subj_path, file_name);

    fprintf('Loading file: %s\n', file_path);
    EEG = pop_loadset('filename', file_name, 'filepath', subj_path);

    % Plot raw 
    eegplot(EEG.data, 'srate', EEG.srate, 'eloc_file', EEG.chanlocs);
else
    fprintf('No *OriginalInterpolated.set file found in %s\n', subj_path);
end

