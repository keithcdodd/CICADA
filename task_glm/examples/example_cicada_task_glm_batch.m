%% Example: run separate first-level CICADA task GLMs over several scans
% Each scan gets its own output. This is not a second-level/group GLM.
% Requires completed CICADA outputs and enough RAM for each 4-D BOLD.

bidsRoot = '/path/to/bids_data';
resultRoot = '/path/to/separate/task_glm_results';

% One row per scan: select Auto or Manual IC classification independently.
scans = {
    'sub-001', 'ses-01', 'taskname_run-01', 'auto';
    'sub-002', 'ses-01', 'taskname_run-01', 'auto';
    'sub-003', 'ses-01', 'taskname_run-01', 'auto'  % change to 'manual' only if completed
};

% Prespecify shared design choices; names must exist in each scan's
% hrf_conditions_raw.csv, and each scan must have matching task timing.
contrasts = struct('name','contrastname', ...
    'conditions',{{'condition1','condition2'}},'weights',[1,-1]); % condition1 - condition2

% Optional dry run: uncomment cfg.dryRun=true below.

for k = 1:size(scans,1)
    cfg = struct();
    cfg.bidsRoot = bidsRoot;
    cfg.subject = scans{k,1};
    cfg.session = scans{k,2};
    cfg.runFolder = scans{k,3};
    cfg.icSelection = scans{k,4};
    cfg.namedContrasts = contrasts;
    cfg.driftModel = 'cosine';
    cfg.highPassHz = 0.008;
    cfg.model = 'AR1';
    cfg.resultRoot = resultRoot;
    % cfg.dryRun = true; % check inputs/design only; no fitting/writes

    R = cicada_task_glm(cfg); %#ok<NASGU>
    fprintf('Scan %d/%d: %s (%s)\n',k,size(scans,1),cfg.subject,R.status);
end

% The output maps are individual-scan first-level effects.
% No across-subject fixed/mixed-effects model is performed by this loop.
