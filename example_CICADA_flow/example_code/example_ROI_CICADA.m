%% Example: independent ROI mean-BOLD time series (no atlas dependency)
% Can be useful for resting-state or task-based connectivity analysis
% Run after an fMRI scan is already processed. ROI extraction is a
% SEPARATE descriptive utility; it does not run ICA, task GLM or Group QC.
% Input BOLD may be subject-level CICADA cleaned or a suitable unsmoothed
% Group-prepared scan. For small anatomical ROIs, prefer non-smoothed data.
% Atlas geometry must already match BOLD: There is no implicit resampling.
%
% 1. EDIT PATHS. ROI files may be in any atlas, not just Brainnetome.
cicadaRepo = '/ABS/PATH/TO/CICADA';
scanDir = '/ABS/PATH/TO/derivatives/cicada/sub-001/ses-01/taskname_run-01';
boldFile = '/ABS/PATH/TO/UNSMOOTHED_CICADA_BOLD.nii.gz';

% Examples: a labeled parcellation (one ROI per positive integer label)
% and optional probabilistic/weighted maps (one ROI per 3-D/4-D volume).
roiFiles = { ...
    '/ABS/PATH/TO/MY_INTEGER_LABELED_ATLAS.nii.gz', ...
    '/ABS/PATH/TO/MY_WEIGHTED_ROI_STACK.nii.gz'};
roiTypes = {'labels', 'weighted'}; % or 'binary' for 3-D/4-D 0/1 masks

resultDir = '/ABS/PATH/TO/NEW_ROI_RESULTS_DIRECTORY';

% Add only the public entry folder; helpers are added/restored internally.
addpath(fullfile(cicadaRepo, 'roi'), '-begin');

% 2. Config: require the scan's canonical quality mask by default.
% This is DISTINCT from CICADA's broad funcmask.nii.gz used for ICA.
cfg = struct();
cfg.boldFile = boldFile;
cfg.roiFiles = roiFiles;
cfg.roiTypes = roiTypes;
cfg.scanDir = scanDir;     % resolves funcmask_constrained.nii.gz
% Alternatively: cfg.qualityMaskFile = '/ABS/PATH/TO/ALIGNED_MASK.nii.gz';
% Omit a quality mask ONLY deliberately: cfg.allowNoQualityMask = true;
cfg.minCoverageFraction = 0; % user-defined, no validated universal cutoff
cfg.minValidVoxels = 1;

% 3. Dry run: validate atlas geometry/types and inspect raw mask coverage.
cfg.dryRun = true;
P = cicada_roi_timeseries(cfg);
disp(P.coverage);

% 4. Deliberately opt in to writing three new files (no overwriting):
% roi_mean_timeseries.csv, roi_coverage.csv, roi_provenance.json.
% The CSV columns roi_001 etc. map to atlas labels in roi_coverage.csv.
% A dry run does not read all BOLD time points; full extraction does.
doExtractAndExport = false;
if doExtractAndExport
    cfg.dryRun = false;
    cfg.writeOutputs = true;
    cfg.outputDir = resultDir;  % parent must exist; directory must NOT exist
    R = cicada_roi_timeseries(cfg);
    fprintf('ROI extraction: %s\n', R.status);
    disp(R.coverage);
end

% NOTE: this utility does not compute ROI contrast beta/variance or
% significance maps. It extracts descriptive, QC-screened BOLD means.
