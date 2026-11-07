function [X, fs, E] = read_eeglab_set(f, varargin)
% Minimal EEGLAB .set/.fdt reader. Handles .set files with top-level fields (saved by pop_saveset) and with an "EEG" struct.
w = whos('-file', f);
if any(strcmp({w.name}, 'EEG'))
    S = load('-mat', f, 'EEG'); S = S.EEG;
else
    S = load('-mat', f, 'data', 'nbchan', 'pnts', 'srate', 'chanlocs');
end
fs = double(S.srate);
if ischar(S.data)
    fid = fopen(fullfile(fileparts(f), S.data), 'r'); X = fread(fid, [double(S.nbchan), double(S.pnts)], 'float32=>double'); fclose(fid);
else
    X = double(S.data);
end
E = struct('nbchan', S.nbchan, 'pnts', S.pnts, 'srate', fs, 'labels', {{}});
if isfield(S, 'chanlocs') && ~isempty(S.chanlocs), E.labels = {S.chanlocs.labels}; end
end
