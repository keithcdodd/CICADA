%% Example: optional CICADA task GLM for ONE completed task scan
% Copy/modify this script. It should be run after Automatic or Manual CICADA.
% COPE and provisional VARCOPE.
% For Group-prepared BOLD, use example_cicada_task_glm_group_matched.m.

cfg = struct();
cfg.bidsRoot = '/path/to/bids_data';
cfg.subject = 'sub-001';
cfg.session = 'ses-01';
cfg.runFolder = 'taskname_run-01';
cfg.icSelection = 'auto';  % 'manual' if Manual CICADA was completed
cfg.inputMode = 'original_joint'; % original BOLD + joint ICA-noise model

cfg.namedContrasts = struct('name','contrastname', ...
    'conditions',{{'condition1','condition2'}},'weights',[1,-1]); % condition1-condition2
cfg.driftModel = 'cosine';  % or 'polynomial'
cfg.highPassHz = 0.008;    % for cosine only; choose before inspecting results
cfg.model = 'AR1';         % or 'AR2' (exploratory covariance models)
cfg.resultRoot = '/path/to/separate/task_glm_results';

% Safe first pass: verify the design without fitting or writing files.
cfg.dryRun = true;
P = cicada_task_glm(cfg);
fprintf('Original-joint task-GLM preflight: %s\n', P.status);

% Deliberate opt-in to a FULL first-level fit and export. This loads the
% entire 4-D image; check memory first. Output must be a NEW location.
doFullFitAndExport = false;
if doFullFitAndExport
    cfg.dryRun = false;
    cfg.runAnalysis = true;
    cfg.writeOutputs = true;
    cfg.maxVoxels = Inf;
    cfg.allowAllVoxels = true;
    cfg.allowLargeRun = true;
    R = cicada_task_glm(cfg);
    fprintf('Original-joint first-level export: %s\n', R.status);
end
