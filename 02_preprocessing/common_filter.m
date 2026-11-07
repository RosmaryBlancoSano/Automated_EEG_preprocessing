function Y = common_filter(X, fs)
% 50 Hz Butterworth band-stop (4th order) + 1 Hz Butterworth high-pass (4th order), zero phase.
Y = zeros(size(X));
wo = 50 / (fs/2); bw = wo / 20;
[bn, an] = butter(4, [wo-bw, wo+bw], 'stop');
[bh, ah] = butter(4, 1/(fs/2), 'high');
for ch = 1:size(X,1)
    Y(ch,:) = filtfilt(bh, ah, filtfilt(bn, an, double(X(ch,:))));
end
end
