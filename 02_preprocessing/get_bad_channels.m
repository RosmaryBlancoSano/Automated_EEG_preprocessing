
function [bad_channels_idx, flat_channels_idx, good_channels_idx] = get_bad_channels(EEG, fpath_out, subj)

    % Create output folders
    subj_out_folder = fullfile(fpath_out, sprintf('sub-%04d', subj));
        if ~exist(subj_out_folder, 'dir')
            mkdir(subj_out_folder);
        end
    figures_folder = fullfile(subj_out_folder, 'Figures');
    if ~exist(figures_folder, 'dir')
        mkdir(figures_folder);
    end

    original_chanlocs = EEG.chanlocs;

    %% --- 1. Bad channel detection via spectral probability
    [EEG_int, ~, measure, ~] = pop_rejchan(EEG, 'elec', 1:EEG.nbchan, ...
        'threshold', [-2.0 2.0], 'measure', 'spec', 'norm', 'on');

    % Plot and save figure
    figure;
    bar(measure);
    xlabel('Channels'); ylabel('Probability');
    title('Probability of being a bad channel per channel');
    saveas(gcf, fullfile(figures_folder, sprintf('bad_chans_sub-%04d.fig', subj)));
    close(gcf);

    % Get bad channels from removedchans
    if isfield(EEG_int.chaninfo, 'removedchans') && ~isempty(EEG_int.chaninfo.removedchans)
        bad_labels = {EEG_int.chaninfo.removedchans.labels};
        bad_channels_idx = find(ismember({original_chanlocs.labels}, bad_labels));
    else
        bad_channels_idx = [];
    end

    %% --- 2. Flat channel detection (low standard deviation)
    channel_stds = std(double(EEG.data), 0, 2);
    flat_threshold = 0.5; 
    flat_channels_idx = find(channel_stds < flat_threshold);

    % Plot
    figure;
    subplot(2, 1, 1);
    bar(channel_stds);
    xlabel('Channels'); ylabel('Standard Deviation (V)');
    title('Standard Deviation per Channel');
    yline(flat_threshold, 'r--', 'Flat threshold');

    subplot(2, 1, 2);
    bar(zscore(channel_stds));
    xlabel('Channels'); ylabel('Z-score of STD');
    title('Z-score of STD (diagnostic)');
    saveas(gcf, fullfile(figures_folder, sprintf('flat_chans_sub-%04d.fig', subj)));
    close(gcf);

    %% --- Combine bad + flat channels
    bad_channels_idx = union(bad_channels_idx, flat_channels_idx);

    %% --- Good channels
    good_channels_idx = setdiff(1:EEG.nbchan, bad_channels_idx);
end
