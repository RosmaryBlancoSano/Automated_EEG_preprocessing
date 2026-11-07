function get_eegdenoisenet(folder)
%GET_EEGDENOISENET  Downloads the public EEGdenoiseNet artefact epochs (EOG_all_epochs.mat, EMG_all_epochs.mat) if they are not in "folder".
%   Source: https://github.com/ncclabsustech/EEGdenoiseNet (folder data/). Needs an internet connection only the first time.
%   If the download is not possible, copy the two files by hand into "folder".
files = {'EOG_all_epochs.mat', 'EMG_all_epochs.mat'};
base = 'https://raw.githubusercontent.com/ncclabsustech/EEGdenoiseNet/master/data/';
if ~isfolder(folder), mkdir(folder); end
for i = 1:numel(files)
    f = fullfile(folder, files{i});
    if ~isfile(f)
        fprintf('Downloading %s from EEGdenoiseNet ...\n', files{i});
        try
            websave(f, [base files{i}]);
        catch ME
            error('get_eegdenoisenet:download', ['Could not download %s (%s).\nDownload it from https://github.com/ncclabsustech/EEGdenoiseNet/tree/master/data ' ...
                  'and copy it into %s'], files{i}, ME.message, folder);
        end
    end
end
end
