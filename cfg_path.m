function p = cfg_path(key)
%CFG_PATH  Folders used by the evaluation scripts (Simulated_groundtruth_evaluation, Real_datasets_evaluation). Edit the PATHS section (or set the environment variables).
%   cfg_path('setup') adds the numbered code folders (01_... to 09_...) of this repository to the MATLAB path.
%   The Stroke and Dementia recordings are patient data and are not part of the repository.
repo = fileparts(mfilename('fullpath'));
% ----------------------------------------------------------------- PATHS (edit)
% Default layout: folders 'data' and 'results' inside the repository (both are ignored by git). If folders with these names exist
% next to the repository folder (the layout of the authors), those are used instead. Each folder can also be set with an
% environment variable, e.g. setenv('EEGDN_RESULTS', 'D:\my_results').
P.results       = getenv_default('EEGDN_RESULTS', pick_dir(repo, 'results'));                                   % csv outputs
P.eegdenoisenet = getenv_default('EEGDN_NET',    fullfile(pick_dir(repo, 'data'), 'SimulatedEEG'));            % public EEGdenoiseNet, downloaded automatically (get_eegdenoisenet.m)
P.auto          = getenv_default('EEGDN_AUTO',   fullfile(pick_dir(repo, 'data'), 'AutomaticArtefactRemoval')); % public: s01 ... s13, each with s01.set
P.dementia      = getenv_default('EEGDN_DEM',    fullfile(pick_dir(repo, 'data'), 'Dementia'));                 % private: ID*/*preproce*.set and *Muse_output.csv
P.stroke        = getenv_default('EEGDN_STROKE', fullfile(pick_dir(repo, 'data'), 'Stroke'));                   % private: sub-XXXX/sub-XXXX_*_denoised.set
P.stroke_input  = getenv_default('EEGDN_STROKE_INPUT', fullfile(pick_dir(repo, 'data'), 'Stroke_GoodQuality'));       % private: GoodQuality recordings (output of 01_quality_check)
P.ids718        = getenv_default('EEGDN_IDS',    fullfile(pick_dir(repo, 'results'), 'ids_718.txt'));           % private: list of the 718 analysed Stroke recordings
% -----------------------------------------------------------------------------
if strcmp(key, 'setup')
    addpath(repo); d = dir(repo);
    for k = 1:numel(d), if d(k).isdir && ~isempty(regexp(d(k).name, '^\d\d_', 'once')), addpath(genpath(fullfile(repo, d(k).name))); end, end
    p = [];
    return
end
p = P.(key);
if strcmp(key, 'results') && ~isfolder(p), mkdir(p); end   % the results folder is created when it is first asked for
end

function v = getenv_default(name, d)
v = getenv(name); if isempty(v), v = d; end
end

function d = pick_dir(repo, name)
sib = fullfile(fileparts(repo), name);
if isfolder(sib), d = sib; else, d = fullfile(repo, name); end
end
