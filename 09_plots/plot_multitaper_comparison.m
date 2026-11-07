
function plot_multitaper_comparison(subj, fpath, originalEEG, denoisedEEG_cell, methodNames, fs, fpath_out)
    % Multitaper settings
    frequency_range = [1 40]; % 1 to 45 Hz
    taper_params = [3, 5]; % time-bandwidth product and number of tapers
    window_params = [4, 1]; % window 4 s, step 1 s
    min_NFFT = 0; % No minimum NFFT
    detrend_opt = 'linear'; % linear detrend
    weighting = 'unity'; % unity weighting of the tapers
    plot_on = true; % plot the spectrograms
    verbose = true; % print log messages

    subj_folder = fullfile(fpath_out, sprintf('Subj%d', subj));
    if ~exist(subj_folder, 'dir')
        mkdir(subj_folder);
    end
    
    % Create the folder if it does not exist
    figures_folder = fullfile(subj_folder, 'figures');
    if ~exist(figures_folder, 'dir')
        mkdir(figures_folder);
    end

    % ---- Multitaper of the original signal ----
    fprintf('Calculating and plotting multitaper for original EEG...\n');
    [original_spect, sfreqs, stimes] = multitaper_spectrogram(originalEEG(:), fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose);

    
    % Save the plot
    saveas(gcf, fullfile(figures_folder, sprintf('Multitaper_Original_Subj%d.png', subj)));

    % ---- Multitaper of every denoised signal ----
    num_methods = numel(denoisedEEG_cell);
    for i = 1:num_methods
        fprintf('Calculating and plotting multitaper for denoised EEG using method: %s...\n', methodNames{i});
        
        % Multitaper spectrogram
        [denoised_spect{i}, sfreqs, stimes] = multitaper_spectrogram(denoisedEEG_cell{i}(:), fs, frequency_range, taper_params, window_params, min_NFFT, detrend_opt, weighting, plot_on, verbose);
        

        % Save the plot
        saveas(gcf, fullfile(figures_folder, sprintf('Multitaper_%s_Subj%d.png', methodNames{i}, subj)));
    
    end
 
end


