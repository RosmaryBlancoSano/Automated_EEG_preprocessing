function data = fake_EEG
% Generate a simulated EEG signal with alpha and beta rhythms
%
% Parameters of the alpha and beta rhythms:
% Ritmo Alpha sugli elettrodi posteriori (O1, O2)
% Ritmo Beta sugli elettrodi fronto-centrali (F3, F4, C3, C4)

% General parameters of the signal
FRAMES = 100;   % frames per trial
TRIALS = 973;   % simulated trials
SRATE = 512;    % sampling rate

% Parameters of the alpha rhythm (posterior channels, e.g. O1, O2)
NEAMP = 75;    % Amplitude for alpha signal (uV)
NEFREQ = 10;    % frequency of the alpha signal (8-13 Hz)
NEPOS = 285;    % position of the alpha peak (in frames), shifted posteriorly

% Parameters of the beta rhythm (fronto-central channels, e.g. F3, F4, C3, C4)
PEAMP = 7;     % Amplitude for beta signal (uV)
PEFREQ = 20;    % frequency of the beta signal (13-30 Hz)
PEPOS = 100;    % position of the beta peak (in frames)

% Jitter of the peaks
TJITTER = 8;    % jitter of the peaks
NOISEAMP = 15; % Amplitude of added noise (uV)

disp('Please wait while generating the data');

% Alpha signal (posterior channels)
alpha_signal = NEAMP * peak(FRAMES, TRIALS, SRATE, NEFREQ, NEPOS, TJITTER);

% Beta signal (fronto-central channels)
beta_signal = PEAMP * peak(FRAMES, TRIALS, SRATE, PEFREQ, PEPOS, TJITTER, 1);

% Load the dipole for the spatial distribution of the signals
Dd = load(fullfile(fileparts(mfilename('fullpath')), 'dipole.mat')); dipole = Dd.dipole;

% Combination of the alpha and beta signals on the appropriate channels
% alpha_signal sui canali posteriori (es. O1, O2)
% beta_signal sui canali fronto-centrali (es. F3, F4, C3, C4)
data = dipole(:, 1) * alpha_signal + dipole(:, 4) * beta_signal;

% Add noise to every channel
for ch = 1:size(data, 1)
    
    data(ch, :) = data(ch, :) + NOISEAMP * noise(FRAMES, TRIALS, SRATE);
end


end

% Generate the peaks (alpha and beta)
function signal = peak(frames, epochs, srate, peakfr, position, tjitter, wave)
% Signal made of a single peak in every trial
% A peak or a sinusoid for every trial

if nargin < 5
   position = frames / 2;
end
if nargin < 6
   tjitter = 0;
end
if nargin < 7
   wave = 0;
end

signal = zeros(1, epochs * frames);
for trial = 1:epochs
   pos = position + round(randn(1) * tjitter);
   for i = 1:frames
       phase = (i - pos) / srate * 2 * pi * peakfr;
       if wave | (phase < pi / 2 & phase > -pi / 2)
           signal((trial - 1) * frames + i) = cos(phase);
       end
   end
end

end

% Generate the noise (simulation of EEG noise)
function signal = noise(frames, epochs, srate)
% Noise with the power spectrum of the human EEG
% It uses a frequency distribution typical of EEG noise

Mm = load(fullfile(fileparts(mfilename('fullpath')), 'meanpower.mat')); meanpower = Mm.meanpower;  
sumsig = 50;  % number of sinusoids that make up every simulated signal

signal = zeros(1, epochs * frames);
for trial = 1:epochs
   freq = 0;
   range = [(trial - 1) * frames + 1 : trial * frames];
   for i = 1:sumsig
      freq = freq + (4 * rand(1));
      freqamp = meanpower(min(ceil(freq), 125)) / meanpower(1);
      phase = rand(1) * 2 * pi;
      signal(range) = signal(range) + sin([1 : frames] / srate * 2 * pi * freq + phase) * freqamp;
   end
end

end


