addpath(fileparts(fileparts(mfilename('fullpath')))); cfg_path('setup');   % folders of the repository (see cfg_path.m)

%% automaticArtifactRemoval
clear; clc;

fpath_out = fullfile(cfg_path('auto'), 'Results');
base_path = cfg_path('auto');

% Define the subject IDs to process
subject_ids = [1:13];

% Loop through each subject
for subj = subject_ids
    % Set file path and file name dynamically
    subj_folder = fullfile(base_path, sprintf('s%02d', subj));
    file_name = sprintf('s%02d.set', subj);
    file_path = fullfile(subj_folder, file_name);
    
    % Load EEG dataset
    EEG_all = pop_loadset('filename', file_name, 'filepath', subj_folder);

    % Extract EEG data and metadata
    data_all = EEG_all.data;
    time_all = EEG_all.times;
    fs_all = EEG_all.srate;
    num_channels_all = size(data_all, 1);

    % Extract channel labels
    Channels = {EEG_all.chanlocs.labels};
    disp('Original Channels:');
    disp(Channels);
    
    %% EEG with 27 channels
    % Exclude EOG channels and update the EEG matrix
    exclude_channels = {'EOG1', 'EOG2', 'EOG3'};
    keep_indices = ~ismember(Channels, exclude_channels);
    data = data_all(keep_indices, :);
    num_channels = size(data, 1);
    Channels = Channels(keep_indices);
    disp(Channels);

    % Update EEG struct to reflect removed channels
    EEG_27.data = data;
    EEG_27.nbchan = num_channels;
    EEG_27.chanlocs = EEG_all.chanlocs(keep_indices); 
    EEG_27.srate = fs_all;
    
    EEG_27_ch = EEG_27.data;   % amplitudes in uV are kept (no z-score): the amplitude-based scores need them
    
    % Save the simulated EEG data inside the subject's folder
    save(fullfile(subj_folder, sprintf('EEG_27_subject%d.mat', subj)), 'EEG_27');
    
    
    %% EEG with 16 channels
    % Exclude channels and update the EEG matrix
    exclude_channels = {'Fz', 'FC5', 'FC1', 'Cz', 'CP5', 'CP1', ...
                    'Pz', 'CP6', 'CP2', 'FC6', 'FC2'};
    keep_indices = ~ismember(Channels, exclude_channels);
    data_temp = EEG_27_ch(keep_indices, :);
    num_channels = size(data_temp, 1);
    Channels_16 = Channels(keep_indices);
    disp(Channels_16);
    
    % Update EEG struct to reflect removed channels
    EEG_16.data = data_temp;
    EEG_16.nbchan = num_channels;
    EEG_16.chanlocs = EEG_27.chanlocs(keep_indices); 
    EEG_16.srate = fs_all;
    
    EEG_16_ch = EEG_16.data;
    
    % Save the simulated EEG data inside the subject's folder
    save(fullfile(subj_folder, sprintf('EEG_16_subject%d.mat', subj)), 'EEG_16');
    
    
    %% EEG with 8 channels
    % Exclude channels and update the EEG matrix
    exclude_channels = {'F3', 'F7', 'P3', 'P7', 'P4', 'P8', 'F4', 'F8'};
    keep_indices = ~ismember(Channels_16, exclude_channels);
    data_temp = EEG_16_ch(keep_indices, :);
    num_channels = size(data_temp, 1);
    Channels_8 = Channels_16(keep_indices);
    disp(Channels_8);
    
    % Update EEG struct to reflect removed channels
    EEG_8.data = data_temp;
    EEG_8.nbchan = num_channels;
    EEG_8.chanlocs = EEG_16.chanlocs(keep_indices); 
    EEG_8.srate = fs_all;
    
    EEG_8_ch = EEG_8.data;
    
    % Save the simulated EEG data inside the subject's folder
    save(fullfile(subj_folder, sprintf('EEG_8_subject%d.mat', subj)), 'EEG_8');
    
    
    %% EEG with 4 channels
    % Exclude channels and update the EEG matrix
    exclude_channels = {'Fp1', 'Fp2', 'F7', 'C3', 'P3', 'P7', 'O1', 'O2', 'P4', 'P8', 'C4', 'F8'};
    keep_indices = ~ismember(Channels_16, exclude_channels);
    data_temp = EEG_16_ch(keep_indices, :);
    num_channels = size(data_temp, 1);
    Channels_4 = Channels_16(keep_indices);
    disp(Channels_4);
    
    % Update EEG struct to reflect removed channels
    EEG_4.data = data_temp;
    EEG_4.nbchan = num_channels;
    EEG_4.chanlocs = EEG_16.chanlocs(keep_indices); 
    EEG_4.srate = fs_all;
    
    EEG_4_ch = EEG_4.data;
    
    % Save the simulated EEG data inside the subject's folder
    save(fullfile(subj_folder, sprintf('EEG_4_subject%d.mat', subj)), 'EEG_4');
    
    
    %% Apply different denoising methods to different chan configs datasets
    
    % Define the file path
    file_names = {sprintf('EEG_4_subject%d.mat', subj), sprintf('EEG_8_subject%d.mat', subj), sprintf('EEG_16_subject%d.mat', subj), sprintf('EEG_27_subject%d.mat', subj)};

    % Define the denoising methods to apply
    denoising_methods = {'W', 'ASR', 'EMD', 'EMD_W', 'ASR_W'};   % methods of the paper (W.m, ASR.m, EMD.m, EMD_W.m, ASR_W.m) % Add more as needed
    % Initialize a structure to store results
    results = struct();
    
    data_loaded = cell(1, length(file_names));
    % Loop su tutti i file di configurazione
    for i = 1:length(file_names)
        % Load the file
        loaded_data = load(fullfile(subj_folder, file_names{i}));  % file of this participant, saved above
        data_loaded{i} = loaded_data;
    end
    
    for j = 1:length(data_loaded)
        loaded_data = data_loaded{j}; % loaded struct
        % EEG field (its name can vary)
        field_names = fieldnames(loaded_data);
        EEG = loaded_data.(field_names{1}); % first field (assumed to be the EEG structure)
     
        % Apply common preprocessing
        fs = EEG.srate;
        filtered_eeg = common_filter(EEG.data, fs);
        EEG.data_filt_orig = filtered_eeg;
        EEG_original_int = EEG; EEG_original_int.data = filtered_eeg;   % input of the method scripts

        % plot
        eegplot(EEG.data_filt_orig, 'srate', EEG.srate, 'eloc_file', EEG.chanlocs);
        close(gcf);
        
        % Loop through denoising methods
        % Loop over the denoising methods
        for k = 1:length(denoising_methods)
            method = denoising_methods{k};

            % Profiling: execution time and resources
            profile on;
            tic;

            % Apply the denoising method
            feval(method);  % call the denoising script by name

            % Cattura i metriche di esecuzione
            elapsedTime = toc;
            profile off;
            p = profile('info');
            cpuTime = sum([p.FunctionTable.TotalTime]);  
            [memUsed, ~] = memory;
            memoryUsed = memUsed.MemUsedMATLAB / 1e6;

            % Check that the denoised field exists and read it
            denoised_field = sprintf('data_%s', method);  % name of the field
            if isfield(EEG, denoised_field)
                data_denoised = EEG.(denoised_field);  % denoised data
            else
                error('Field "%s" does not exist in the EEG structure. Check the denoising function.', denoised_field);
            end

            % Store the results
            fieldname = sprintf('Subj%d_%s_%s', subj, field_names{1}, strrep(method, '_', '')); 
        
            results.(fieldname).data_denoised = data_denoised;
            results.(fieldname).elapsedTime = elapsedTime;
            results.(fieldname).cpuTime = cpuTime;
            results.(fieldname).memoryUsed = memoryUsed;

            % Display the results
            fprintf('Subject: %d, Configuration: %s, Method: %s\nElapsed time: %.2f s\nCPU time: %.2f s\nMemory used: %.2f MB\n\n', ...
                subj, field_names{1}, method, elapsedTime, cpuTime, memoryUsed);

            % Create the folder of the participant if it does not exist
            subj_folder = fullfile(fpath_out, sprintf('Subj%d', subj));
            if ~exist(subj_folder, 'dir')
                mkdir(subj_folder);
            end


            h = eegplot(EEG.data_filt_orig, 'data2', results.(fieldname).data_denoised, 'events', EEG_all.event, 'title', method);
            close(gcf);

            % Save the plot as .fig in the folder of the participant
        end

        % Cell array denoisedEEG_cell with the denoised signals of all methods denoised
        denoisedEEG_cell = cell(1, length(denoising_methods));
        for i = 1:length(denoising_methods)
            method = denoising_methods{i};
            denoised_field = sprintf('data_%s', method); % name of the field
            if isfield(EEG, denoised_field)
                denoisedEEG_cell{i} = EEG.(denoised_field);
            else
                error('Field "%s" does not exist in the EEG structure. Check the denoising function.', denoised_field);
            end
        end

        % Call compareDenoisingMethods for post-processing
        originalEEG = EEG.data_filt_orig; 
        methodNames = denoising_methods;
        d_events = {EEG_all.event.type};
        lambda = 0.5; 
        it = length(denoising_methods); % Number of iterations/methods
        fpath = file_path; % Define file path for function input
        channelConfig = field_names{1};
        
        compareDenoisingMethods(subj, originalEEG, denoisedEEG_cell, methodNames, lambda, fs, fpath_out, channelConfig);  
        

        
        saveDenoisingPerformance(results, fpath_out, channelConfig);

    end
end



