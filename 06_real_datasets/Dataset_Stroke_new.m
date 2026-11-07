addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');   % folders of the repository (see cfg_path.m)
%% ---------------------------------------------------------------------
%% Dataset: Stroke
%% ---------------------------------------------------------------------

%% Set Paths
clear; clc;
% EEGLAB (https://sccn.ucsd.edu/eeglab, 2022.1 or later) must be on the MATLAB path and started (eeglab nogui) before running

fpath_out = cfg_path('stroke');
base_path = cfg_path('stroke_input');

% Subjects Folder list which start with 'sub-'
subj_folders = dir(fullfile(base_path, 'sub-*'));
subj_folders = subj_folders([subj_folders.isdir]); 

% Loop through subject folders
for i = 1:length(subj_folders)
    subj_folder = subj_folders(i).name;
    subj_path = fullfile(base_path, subj_folder);

    % Extract subject number from folder name
    subj = str2double(regexp(subj_folder, '\d+', 'match', 'once')); 
    
    % Look for the *.set file
    eeg_file = dir(fullfile(subj_path, '*.set'));
    
    if ~isempty(eeg_file)
        file_name = eeg_file(1).name;
        file_path = fullfile(subj_path, file_name);

        fprintf('Loading file: %s\n', file_path);
        EEG = pop_loadset('filename', file_name, 'filepath', subj_path);

        num_samples_per_second = EEG.srate;
        total_samples = EEG.pnts;
        samples_to_remove = num_samples_per_second;  % 1 second in samples
        start_sample = samples_to_remove;  % first second to remove
        end_sample = total_samples - samples_to_remove;  % last second to remove
        
        % remove the first and the last second of the recording
        EEG.data = EEG.data(:, start_sample:end_sample);

        disp(EEG.pnts)      % number of samples
        disp(EEG.srate)     % sampling rate
        disp(EEG.xmax)      % duration in seconds
        size(EEG.data)      % dimensions [nchan, npnts]

        EEG = pop_chanedit(EEG, 'lookup', fullfile(fileparts(which('eeglab')), 'sample_locs', 'Standard-10-20-Cap8.ced'));
        
        %% Extract channel labels
        Channels = {EEG.chanlocs.labels};
        disp('Original Channels:');
        disp(Channels);

    end

%% 1. Filtering
    % The Stroke recordings are already filtered (FIR notch at 25, 50 and 100 Hz, 1 Hz high-pass, 35 Hz low-pass) and cleaned of
    % high-amplitude segments by 01_quality_check/DataQualityCheck.py; no further filtering is applied here
    % (common_filter is used for the semi-simulated and AutomaticArtefactRemoval datasets).

    % Save and Plot filtered EEG

    subj_out_folder = fullfile(fpath_out, sprintf('sub-%04d', subj));
        if ~exist(subj_out_folder, 'dir')
            mkdir(subj_out_folder);
        end

    figures_folder = fullfile(subj_out_folder, 'Figures');  % folder for the figures
        if ~exist(figures_folder, 'dir')
            mkdir(figures_folder);
        end

        pop_saveset(EEG, ...
        'filename', sprintf('sub-%04d_Original.set', subj), ...
        'filepath', subj_out_folder);

        eegplot(EEG.data, 'srate', EEG.srate, 'title', 'Original', 'eloc_file', EEG.chanlocs);
        saveas(gcf, fullfile(figures_folder, sprintf('eegplot_sub-%04d_Original.fig', subj)));
        close(gcf);
        
    %% 2. Detect and interpolate bad channels:

        EEG_copy = EEG;

        [bad_channels_idx, ~, ~] = get_bad_channels(EEG_copy, base_path, subj);
    
        fprintf('Subject %04d: %d bad channels detected\n', subj, length(bad_channels_idx));
    
        if ~isempty(bad_channels_idx)
            EEG_original_int = eeg_interp(EEG_copy, bad_channels_idx, 'spherical');
        else
            EEG_original_int = EEG_copy; % no interpolation needed
        end


        %% Save Origina lnterpolated EEG Plot in subjec folder
        eegplot(EEG_original_int.data, 'srate', EEG.srate, 'title', 'Original Interp', 'eloc_file', EEG.chanlocs);
        saveas(gcf, fullfile(figures_folder, sprintf('eegplot_sub-%04d_Original_Interp.fig', subj)));
        close(gcf);

        pop_saveset(EEG_original_int, ...
        'filename', sprintf('sub-%04d_OriginalInterpolated.set', subj), ...
        'filepath', subj_out_folder);
        
    %% 3. Apply different denoising methods

    denoising_methods = {'W', 'ASR', 'EMD',  'EMD_W', 'ASR_W'};  % , 'EMD_WT_NEW'
    results = struct();  

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
        output_file = fullfile(subj_out_folder, sprintf('sub-%04d_%s_denoised.set', subj, method));
        EEG.data = data_denoised;  % Assign the denoised data to the EEG structure
        EEG.setname = sprintf('sub-%04d_%s_denoised.set', subj, method);  % Set a name for the set
        pop_saveset(EEG, 'filename', output_file);  % Save EEG data as .set file

        % if method ASR, apply cross-corr to correcting for lag
        if strcmp(method, 'ASR')

            [xc, lags] = xcorr(EEG_original_int.data(1,:), results.(fieldname).data_denoised(1,:), 'coeff'); 
            [~, maxIdx] = max(xc);
            b = lags(maxIdx);
            data_denoised_shifted = circshift(results.(fieldname).data_denoised, b, 2);

            % Plot original vs denoised data (shifted)
            eegplot(EEG_original_int.data, 'data2', data_denoised_shifted, 'srate', EEG.srate, 'title', sprintf('%s (ASR)', method), 'eloc_file', EEG.chanlocs);
            fig1 = gcf;  
            saveas(fig1, fullfile(figures_folder, ...
            sprintf('eegplot_sub-%04d_%s_superimposed.fig', subj, method)));
            close(fig1);

            eegplot(data_denoised_shifted, 'srate', EEG.srate, 'title', sprintf('%s (ASR)', method), 'eloc_file', EEG.chanlocs);
            fig2 = gcf;
            saveas(fig2, fullfile(figures_folder, ...
                   sprintf('eegplot_sub-%04d_%s_denoised.fig', subj, method)));
            close(fig2);

        else
            % Plot original vs denoised data
            eegplot(EEG_original_int.data, 'data2', results.(fieldname).data_denoised, 'srate', EEG.srate, 'title', method, 'eloc_file', EEG.chanlocs);
            fig1 = gcf;  % the figure just created
            saveas(fig1, fullfile(figures_folder, ...
            sprintf('eegplot_sub-%04d_%s_superimposed.fig', subj, method)));
            close(fig1);

            eegplot(results.(fieldname).data_denoised, 'srate', EEG.srate, 'title', method, 'eloc_file', EEG.chanlocs);
            fig2 = gcf;
            saveas(fig2, fullfile(figures_folder, ...
                   sprintf('eegplot_sub-%04d_%s_denoised.fig', subj, method)));
            close(fig2);
        end
    end

    % Store denoised EEG data in a cell array
    denoisedEEG_cell = cell(1, length(denoising_methods));
    for i = 1:length(denoising_methods)
        method = denoising_methods{i};
        denoised_field = sprintf('data_%s', method); 
        method_field = matlab.lang.makeValidName(sprintf('Subj%d_%s', subj, strrep(method, '_', '')));
        if strcmp(method, 'ASR')
            [xc, lags] = xcorr(EEG_original_int.data(1,:), results.(method_field).data_denoised(1,:), 'coeff');
            [~, maxIdx] = max(xc);
            b = lags(maxIdx);
            denoisedEEG_cell{i} = circshift(results.(method_field).data_denoised, b, 2);
        elseif isfield(EEG, denoised_field)
            denoisedEEG_cell{i} = EEG.(denoised_field);
        else
            warning('Field "%s" not found. Skipping...', denoised_field);
        end
    end

    %% 4. Comparison and evaluation 
    % Call compareDenoisingMethods for post-processing
    originalEEG = EEG_original_int.data; 
    methodNames = denoising_methods;
    lambda = 0.5; 
    it = length(denoising_methods); % Number of iterations/methods
    
    compareDenoisingMethods(subj, originalEEG, denoisedEEG_cell, methodNames, lambda, fs, fpath_out, '8ch');

    % Plot and Save time-frequency results
    plot_multitaper_comparison(subj, [], originalEEG, denoisedEEG_cell, methodNames, fs, fpath_out);
    
    % Save Denoised Performance
    saveDenoisingPerformance(results, fpath_out, '8ch');
        
end


