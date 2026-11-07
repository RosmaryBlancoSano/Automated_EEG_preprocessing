addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');   % folders of the repository (see cfg_path.m)
%% ---------------------------------------------------------------------
%% Dataset: Dementia (Muse, 4 channels, 256 Hz)
%% ---------------------------------------------------------------------
% Input : data/Dementia/ID*/*preproce*.set   (TP9, AF7, AF8, TP10 and the marker channel 'stim'; 24 participants, ID0013 excluded)
% Output: data/Dementia/EEG_denoised/IDxxxx/ (original and denoised .set files, figures, performance) and Subjx/ (comparison with the original, performance)
% Same steps as Dataset_Stroke_new.m, except that the recordings are used in full (no seconds removed) and that no channel is
% interpolated: with 4 electrodes spatial interpolation was not considered reliable.

%% Set Paths
clear; clc;
% EEGLAB (https://sccn.ucsd.edu/eeglab, 2022.1 or later) must be on the MATLAB path and started (eeglab nogui) before running

base_path = cfg_path('dementia');
fpath_out = fullfile(base_path, 'EEG_denoised');

% Subjects folder list which start with 'ID'
subj_folders = dir(fullfile(base_path, 'ID*'));
subj_folders = subj_folders([subj_folders.isdir]);
subj_folders = subj_folders(~strcmp({subj_folders.name}, 'ID0013'));   % 24 participants as in the paper

% Loop through subject folders
for i = 1:length(subj_folders)
    subj_folder = subj_folders(i).name;
    subj_path = fullfile(base_path, subj_folder);

    % Extract subject number from folder name
    subj = str2double(regexp(subj_folder, '\d+', 'match', 'once'));

    % Look for the *.set file
    eeg_file = dir(fullfile(subj_path, '*preproce*.set'));
    if isempty(eeg_file)
        warning('No *.set file in %s. Skipping...', subj_path);
        continue;
    end
    fprintf('Loading file: %s\n', fullfile(subj_path, eeg_file(1).name));
    EEG = pop_loadset('filename', eeg_file(1).name, 'filepath', subj_path);
    fs = EEG.srate;

    %% Keep the 4 EEG channels; the marker channel is stored in EEG.stim
    stim_index = find(strcmpi({EEG.chanlocs.labels}, 'stim'), 1);
    eeg_indices = setdiff(1:EEG.nbchan, stim_index);
    EEG.stim = EEG.data(stim_index, :);
    EEG.data = EEG.data(eeg_indices, :);
    EEG.chanlocs = EEG.chanlocs(eeg_indices);
    EEG.nbchan = length(eeg_indices);
    disp('Channels:');
    disp({EEG.chanlocs.labels});

    %% 1. Filtering
    % The recordings are already filtered and cleaned of high-amplitude segments; no further filtering is applied here.

    subj_out_folder = fullfile(fpath_out, sprintf('ID%04d', subj));
    if ~exist(subj_out_folder, 'dir')
        mkdir(subj_out_folder);
    end
    figures_folder = fullfile(subj_out_folder, 'Figures');
    if ~exist(figures_folder, 'dir')
        mkdir(figures_folder);
    end

    % Save and Plot Original EEG
    EEG.setname = sprintf('ID%04d_Original.set', subj);
    pop_saveset(EEG, ...
        'filename', sprintf('ID%04d_Original.set', subj), ...
        'filepath', subj_out_folder);

    eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Original', 'eloc_file', EEG.chanlocs);
    saveas(gcf, fullfile(figures_folder, sprintf('eegplot_ID%04d_Original.fig', subj)));
    close(gcf);

    %% 2. No channel interpolation
    % The denoising scripts read EEG_original_int.
    EEG_original_int = EEG;

    %% 3. Apply different denoising methods

    denoising_methods = {'W', 'ASR', 'EMD', 'EMD_W', 'ASR_W'};
    results = struct();
    denoisedEEG_cell = cell(1, length(denoising_methods));

    % Loop through denoising methods
    for k = 1:length(denoising_methods)
        method = denoising_methods{k};

        % Check if function exists before calling it
        if exist(method, 'file') ~= 2
            warning('Function "%s" not found. Skipping...', method);
            continue;
        end

        % Capture performance metrics
        profile on;
        tic;
        feval(method);
        elapsedTime = toc;
        profile off;
        p = profile('info');
        cpuTime = sum([p.FunctionTable.TotalTime]);
        [memUsed, ~] = memory;
        memoryUsed = memUsed.MemUsedMATLAB / 1e6;

        % Check if denoised field exists in EEG structure
        denoised_field = sprintf('data_%s', method);
        if isfield(EEG, denoised_field)
            data_denoised = EEG.(denoised_field);
        else
            error('Field "%s" does not exist in EEG structure. Check denoising function.', denoised_field);
        end
        denoisedEEG_cell{k} = data_denoised;

        %% Save performance metrics results
        fieldname = matlab.lang.makeValidName(sprintf('Subj%d_%s', subj, strrep(method, '_', '')));
        disp(fieldname)

        results.(fieldname).data_denoised = data_denoised;
        results.(fieldname).elapsedTime = elapsedTime;
        results.(fieldname).cpuTime = cpuTime;
        results.(fieldname).memoryUsed = memoryUsed;

        % Display performance
        fprintf('Subject: %d, Method: %s\nElapsed Time: %.2f s\ncpuTime: %.2f s\nMemory Used: %.2f MB\n\n', ...
            subj, method, elapsedTime, cpuTime, memoryUsed);

        %% Save denoised data as EEGLAB .set file
        output_file = fullfile(subj_out_folder, sprintf('ID%04d_%s_denoised.set', subj, method));
        EEG.data = data_denoised;  % Assign the denoised data to the EEG structure
        EEG.setname = sprintf('ID%04d_%s_denoised.set', subj, method);
        pop_saveset(EEG, 'filename', output_file);

        % Plot original vs denoised data
        eegplot(EEG_original_int.data, 'data2', data_denoised, 'srate', EEG.srate, 'title', method, 'eloc_file', EEG.chanlocs);
        fig1 = gcf;
        saveas(fig1, fullfile(figures_folder, sprintf('eegplot_ID%04d_%s_superimposed.fig', subj, method)));
        close(fig1);

        eegplot(data_denoised, 'srate', EEG.srate, 'title', method, 'eloc_file', EEG.chanlocs);
        fig2 = gcf;
        saveas(fig2, fullfile(figures_folder, sprintf('eegplot_ID%04d_%s_denoised.fig', subj, method)));
        close(fig2);
    end

    %% 4. Comparison and evaluation
    originalEEG = EEG_original_int.data;
    methodNames = denoising_methods;
    lambda = 0.5;

    compareDenoisingMethods(subj, originalEEG, denoisedEEG_cell, methodNames, lambda, fs, fpath_out, '4ch');

    % Plot and Save time-frequency results
    plot_multitaper_comparison(subj, [], originalEEG, denoisedEEG_cell, methodNames, fs, fpath_out);

    % Save Denoised Performance
    saveDenoisingPerformance(results, fpath_out, '4ch');

end
