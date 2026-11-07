
%% SNR_event (Signal-to-Noise Ratio) measures the ratio between the signal power (the useful part) and the noise power (the unwanted part).
% The higher the SNR, the better the denoised signal is compared to the noise. 
% A lower SNR indicates that the noise is stronger than the useful signal. 
% SNR is expressed in decibels (dB), and generally, higher values are better. For example, an SNR of 20 dB means that the signal is 100 times more powerful than the noise, while an SNR of 0 dB means the signal is equal to the noise.
%% PeakSNR_event (Peak Signal-to-Noise Ratio) is an improved version of SNR that considers the maximum value of the signal (not just the average) relative to the noise.
% Higher values of PeakSNR indicate better signal quality. 
% A lower PeakSNR could suggest that the denoised signal is not very good compared to the peaks of the original signal.
%% rrmse_temporal_values: The Root Relative Mean Square Error (RRMSE) in the temporal domain measures the difference between the original and denoised EEG signals over time. A lower RRMSE value indicates that the denoised signal is more similar to the original signal, meaning better denoising performance.
%% rrmse_spectral_values: The Root Relative Mean Square Error (RRMSE) in the spectral domain measures the difference between the original and denoised signals in terms of their frequency components. A lower RRMSE spectral value indicates better preservation of the frequency characteristics after denoising.
%% cc_values: The Pearson Correlation Coefficient (CC) measures the linear correlation between the original and denoised EEG signals. A value of 1 indicates perfect positive correlation, 0 means no correlation, and -1 indicates perfect negative correlation. Higher values represent better alignment between the original and denoised signals.
%% deltaRatio_values, thetaRatio_values, alphaRatio_values, betaRatio_values, gammaRatio_values: This is the ratio of the band power in the denoised signal to the band power in the original signal. Higher values indicate better preservation of band power in the denoised signal.
% In general: higher ratios for each frequency band indicate that the denoising method preserves the respective band power effectively. 
% Lower RRMSE and higher correlation (CC) suggest better performance in maintaining the original signal characteristics.

function resultsComparisonMethods = compareDenoisingMethods(subj, originalEEG, denoisedEEG_cell, methodNames, lambda, fs, fpath_out, channelConfig)
    if nargin < 7
        lambda = 1; % Default lambda value if not provided
    end

    numMethods = length(denoisedEEG_cell);
    resultsComparisonMethods = struct();

    % Windowing parameters
    windowDuration = 2; % window duration in seconds
    overlap = 0.5; % 50% overlap
    windowSize = windowDuration * fs; % samples per window
    stepSize = floor(windowSize * (1 - overlap)); % window step

    % make sure that denoisedEEG_cell is a cell array
    if ~iscell(denoisedEEG_cell)
        denoisedEEG_cell = {denoisedEEG_cell}; % convert to cell array if needed
    end

    for i = 1:numMethods
        denoisedEEG = denoisedEEG_cell{i};
        methodName = methodNames{i};
        resultsComparisonMethods(i).Method = methodName;
        event_results = struct(); % results of every method

        % Analyse the whole signal, divided into windows
        snr_values = [];
        peakSNR_values = [];
        rrmse_temporal_values = [];
        rrmse_spectral_values = [];
        cc_values = [];
        deltaRatio_values = [];
        thetaRatio_values = [];
        alphaRatio_values = [];
        betaRatio_values = [];
        gammaRatio_values = [];

        % Divide the signal into windows
        numWindows = floor((size(originalEEG, 2) - windowSize) / stepSize) + 1;

        for startIdx = 1:stepSize:(size(originalEEG, 2) - windowSize + 1)
            endIdx = startIdx + windowSize - 1;
            if endIdx > size(originalEEG, 2)
                endIdx = size(originalEEG, 2);
            end

            % EEG segments
            original_segment = originalEEG(:, startIdx:endIdx);
            denoised_segment = denoisedEEG(:, startIdx:endIdx);
            noise_segment = original_segment - denoised_segment;

            % Quantitative measures
            [~, snr, rrmseTemporal, rrmseSpectral, cc] = signal_analysis(original_segment, denoised_segment, noise_segment, lambda, fs);

            % SNR and PeakSNR with calcSNR_PSNR
            [~, PeakSNR] = calcSNR_PSNR(original_segment, denoised_segment, 1); % preEEG - postEEG

            % PSD bands
            [originalBands, ~] = calculatePSDBands(original_segment, fs);
            [denoisedBands, ~] = calculatePSDBands(denoised_segment, fs);

            % Add the results of the measures
            snr_values(end + 1) = snr;
            peakSNR_values(end + 1) = PeakSNR;
            rrmse_temporal_values(end + 1) = rrmseTemporal;
            rrmse_spectral_values(end + 1) = rrmseSpectral;
            cc_values(end + 1) = cc;
            deltaRatio_values(end + 1) = denoisedBands.delta / (originalBands.delta + eps);
            thetaRatio_values(end + 1) = denoisedBands.theta / (originalBands.theta + eps);
            alphaRatio_values(end + 1) = denoisedBands.alpha / (originalBands.alpha + eps);
            betaRatio_values(end + 1) = denoisedBands.beta / (originalBands.beta + eps);
            gammaRatio_values(end + 1) = denoisedBands.gamma / (originalBands.gamma + eps);
        end

        % Average of the results of every event type
        event_results.('all') = struct(...
            'SNR', mean(snr_values), ...
            'PeakSNR', mean(peakSNR_values), ...
            'RRMSE_Temporal', mean(rrmse_temporal_values), ...
            'RRMSE_Spectral', mean(rrmse_spectral_values), ...
            'PearsonCC', mean(cc_values), ...
            'DeltaRatio', mean(deltaRatio_values), ...
            'ThetaRatio', mean(thetaRatio_values), ...
            'AlphaRatio', mean(alphaRatio_values), ...
            'BetaRatio', mean(betaRatio_values), ...
            'GammaRatio', mean(gammaRatio_values));

        % Store the results of every method
        resultsComparisonMethods(i).Events = event_results;
    end

    % Save the results in a CSV file
    saveResultsToCSV(resultsComparisonMethods, fpath_out, subj, channelConfig);
    plotResults(resultsComparisonMethods, fpath_out, subj, channelConfig);
end

function saveResultsToCSV(resultsComparisonMethods, fpath_out, subj, channelConfig)
    % Create the folder of the participant if it does not exist
    subj_folder = fullfile(fpath_out, sprintf('Subj%d', subj));
    if ~exist(subj_folder, 'dir')
        mkdir(subj_folder);
    end

    % Header of the CSV file, including the channel configuration
    header = {'Method', 'Event', 'SNR', 'PeakSNR', 'RRMSE_Temporal', 'RRMSE_Spectral', ...
              'PearsonCC', 'DeltaRatio', 'ThetaRatio', 'AlphaRatio', 'BetaRatio', 'GammaRatio', 'ChannelConfig'};
    
    % Cell array for the data
    dataTable = cell(0, length(header));

    % Number of methods and event types
    numMethods = length(resultsComparisonMethods);
    eventTypes = fields(resultsComparisonMethods(1).Events);

    for i = 1:numMethods
        methodName = resultsComparisonMethods(i).Method;
        for j = 1:length(eventTypes)
            eventType = eventTypes{j};
            eventData = resultsComparisonMethods(i).Events.(eventType);
            
            % Row of data for the CSV, including the channel configuration
            rowData = {methodName, eventType, eventData.SNR, eventData.PeakSNR, eventData.RRMSE_Temporal, ...
                       eventData.RRMSE_Spectral, eventData.PearsonCC, eventData.DeltaRatio, ...
                       eventData.ThetaRatio, eventData.AlphaRatio, eventData.BetaRatio, eventData.GammaRatio, channelConfig};
            
            % Add the row to the table
            dataTable(end+1, :) = rowData;
        end
    end

    % Convert the cell array to a table
    T = cell2table(dataTable, 'VariableNames', header);

    % Save the CSV file in the folder of the participant
    csv_name = sprintf('resultsComparisonMethods_Subj%d_%s.csv', subj, channelConfig);
    writetable(T, fullfile(subj_folder, csv_name));

    fprintf('Results saved for Subj%d in %s\n', subj, csv_name);
end

function plotResults(resultsComparisonMethods, fpath_out, subj, channelConfig)
    % Create the folder of the participant if it does not exist
    subj_folder = fullfile(fpath_out, sprintf('Subj%d', subj));
    if ~exist(subj_folder, 'dir')
        mkdir(subj_folder);
    end

    % Subfolder "figures" for the images
    figures_folder = fullfile(subj_folder, 'figures');
    if ~exist(figures_folder, 'dir')
        mkdir(figures_folder);
    end

    % Number of methods
    numMethods = length(resultsComparisonMethods);
    metricNames = {'SNR', 'PeakSNR', 'RRMSE_Temporal', 'RRMSE_Spectral', 'PearsonCC', 'DeltaRatio', 'ThetaRatio', 'AlphaRatio', 'BetaRatio', 'GammaRatio'};

    figure;
    sgtitle(sprintf('Comparison of Denoising Methods - Subj%d - %s', subj, channelConfig));

    for m = 1:length(metricNames)
        metricName = metricNames{m};
        metricValues = zeros(1, numMethods);
        for i = 1:numMethods
            metricValues(i) = resultsComparisonMethods(i).Events.all.(metricName);
        end

        % Normalise the metric values
        if max(metricValues) ~= 0
            normalizedValues = metricValues / max(metricValues);
        else
            normalizedValues = zeros(size(metricValues));
        end

        subplot(4, 3, m);
        colors = lines(numMethods); % default colour palette
        b = bar(1:numMethods, metricValues, 'FaceColor', 'flat'); 

        % a different colour for every bar
        for k = 1:numMethods
            b.CData(k, :) = colors(k, :);
        end

        xticks(1:numMethods);
        xticklabels({resultsComparisonMethods(:).Method});
        xtickangle(45);
        ylabel('Norm value');
        title(metricName);
    end

    % Save the figure in the folder "figures" of the participant
    fig_name = sprintf('Comparison_Subj%d_%s.fig', subj, channelConfig);
    savefig(gcf, fullfile(figures_folder, fig_name));

    fprintf('Figures saved for Subj%d in %s\n', subj, figures_folder);
end

function [rms_val, snr, rrmse_temporal, rrmse_spectral, cc] = signal_analysis(original_segment, denoised_segment, noise_segment, lambda, fs)
    % RMS value
    rms_val = sqrt(mean(original_segment.^2, 2));  

    % SNR
    rms_noise = sqrt(mean((lambda * noise_segment).^2, 2));  
    snr = 10 * log10(mean(rms_val) ./ (mean(rms_noise) + eps));  % mean of the SNR values
    
    % RRMSE in the time domain
    rrmse_temporal = mean(sqrt(mean((denoised_segment - original_segment).^2, 2)) ./ rms_val); 

    % RRMSE in the frequency domain
    [pxx_original, ~] = pwelch(original_segment', [], [], [], fs);
    [pxx_denoised, ~] = pwelch(denoised_segment', [], [], [], fs);
    rrmse_spectral = mean(sqrt(mean((pxx_denoised - pxx_original).^2, 2)) ./ sqrt(mean(pxx_original.^2, 2)));  

    % Pearson correlation coefficient
    
    % Cross-correlation
    cc_values = zeros(1, size(original_segment, 1)); % correlation values
    lag_values = zeros(1, size(original_segment, 1)); % lag values

    for i = 1:size(original_segment, 1)
        [xc, lags] = xcorr(original_segment(i, :), denoised_segment(i, :), 'coeff');
        [~, maxIdx] = max(xc);
        cc_values(i) = xc(maxIdx); % maximum of the cross-correlation
        lag_values(i) = lags(maxIdx); % lag of the maximum
    end

    cc = mean(cc_values); % mean cross-correlation (between -1 and 1)
    lag = mean(lag_values)/fs; % mean lag in s
    
end

function [bands, freqs] = calculatePSDBands(eeg, fs)
    % Compute Power Spectral Density (PSD) using Welch's method
    [psd, freqs] = pwelch(eeg', [], [], [], fs);

    % Compute band powers for different frequency ranges
    bands.delta = bandpower(psd, fs, [1 4]);    % Delta band: 1-4 Hz
    bands.theta = bandpower(psd, fs, [4 8]);    % Theta band: 4-8 Hz
    bands.alpha = bandpower(psd, fs, [8 13]);   % Alpha band: 8-13 Hz
    bands.beta  = bandpower(psd, fs, [13 30]);  % Beta band: 13-30 Hz
    bands.gamma = bandpower(psd, fs, [30 45]);  % Gamma band: 30-45 Hz
end

    
    
