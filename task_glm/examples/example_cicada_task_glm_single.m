%% Example: optional CICADA task GLM for ONE completed task scan
% Copy/modify this script. It runs after Automatic or Manual CICADA.
% Exploratory COPE and provisional VARCOPE only; no t/z/p significance.

cfg = struct();
cfg.bidsRoot = '/path/to/bids_data';
cfg.subject = 'sub-001';
cfg.session = 'ses-01';
cfg.runFolder = 'taskname_run-01';
cfg.icSelection = 'auto';  % 'manual' if Manual CICADA was completed

cfg.namedContrasts = struct('name','contrastname', ...
    'conditions',{{'condition1','condition2'}},'weights',[1,-1]); % condition1-condition2
cfg.driftModel = 'cosine';  % or 'polynomial'
cfg.highPassHz = 0.008;    % for cosine only; choose before inspecting results
cfg.model = 'AR1';         % or 'AR2' (exploratory covariance models)
cfg.resultRoot = '/path/to/separate/task_glm_results';

% One call: validate, fit entire scan, and save results to a NEW folder.
R = cicada_task_glm(cfg);

% Optional inspection without fitting or creating directories:
% cfg.dryRun = true;
% Rpreview = cicada_task_glm(cfg);
