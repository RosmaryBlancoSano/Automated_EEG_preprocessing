function start_pool(n_workers)
%START_POOL  Starts a pool of n_workers local workers if Parallel Computing Toolbox is available; otherwise (or with n_workers = 0)
%   the parfor loops of the evaluation scripts (Simulated_groundtruth_evaluation, Real_datasets_evaluation) run serially on the client.
if n_workers <= 0, return; end
try
    if isempty(gcp('nocreate')), parpool('Processes', n_workers); end
catch
    warning('start_pool:noPool', 'No parallel pool available: running serially.');
end
end
