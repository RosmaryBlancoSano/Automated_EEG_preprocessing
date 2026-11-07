
function denoisedEEG = wavelet_denoising_default(eegData, srate, paradigmERP, softThresh)
    % WAVELET_DENOISING_DEFAULT - Perform wavelet denoising on EEG data.
    % This function applies wavelet-based denoising with configurable
    % parameters.
    %
    % Parameters:
    %   eegData     - EEG data (channels x samples).
    %   srate       - Sampling rate of the EEG signal.
    %   paradigmERP - Boolean, set to true for ERP paradigms.
    %   softThresh  - Boolean, set to true for soft thresholding.
    %
    % Returns:
    %   denoisedEEG - Denoised EEG data (channels x samples).

    % Select wavelet family and decomposition level based on sampling rate and paradigm
    if paradigmERP
        wavFam = 'coif4';
        if srate > 500
            wavLvl = 10;
        elseif srate > 250 && srate <= 500
            wavLvl = 9;
        else
            wavLvl = 8;
        end
    else
        wavFam = 'coif4'; %'bior4.4'; % 'db4'; %
        if srate > 500
            wavLvl = 8;
        elseif srate > 250 && srate <= 500
            wavLvl = 7;
        else
            wavLvl = 6;
        end
    end

    % Set threshold rule
    if paradigmERP && softThresh
        ThresholdRule = 'Soft';
    else
        ThresholdRule = 'Hard';
    end

    % Ensure EEG data is double
    if ~isa(eegData, 'double')
        eegData = double(eegData);
    end

    % Initialize matrix for denoised data
    [numChannels, ~] = size(eegData);
    denoisedEEG = zeros(size(eegData));

    % Process each channel independently
    for ch = 1:numChannels
        % Reshape data if needed (flatten for processing)
        signal = eegData(ch, :);

        % Perform wavelet denoising
        artifacts = wdenoise(signal, wavLvl, 'Wavelet', wavFam, 'DenoisingMethod', 'Bayes', ...
                             'ThresholdRule', ThresholdRule, 'NoiseEstimate', 'LevelDependent');

        % Subtract the artifacts from the original signal
        denoisedEEG(ch, :) = signal - artifacts;
    end

end


