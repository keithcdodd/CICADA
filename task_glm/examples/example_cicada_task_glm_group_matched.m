%% Example: experimental task GLM after verified matched Group CICADA
% This analyzes ONE previously completed task-fMRI scan. It does NOT run
% Auto/Manual/Group CICADA, MELODIC, or second-level inference.
%
% REQUIREMENTS FOR THIS SPECIFIC MODE:
%   - The Group run used temporal_cleanup_mode = 'matched'.
%   - Its preparation settings were detrend_degree = 2,
%     smoothing_kernel = 0 mm and voxelwise_scale_flag = 0.
%   - It is built to handle high pass filtering (take it into account)
%   - Its producer-time *_matched_recipe.json is present and valid.
%   - The source scan has hrf_conditions_raw.csv and verified task timing.
% Older/unattested Group NIfTIs cannot be used as if they were verified.
% The resulting COPE maps and plug-in VARCOPE maps are EXPERIMENTAL;
% no calibrated t, z, p or second-level inference is produced here.
%
% 1. EDIT THESE PATHS / RUN IDENTIFIERS FOR YOUR DATA.
cicadaRepo = '/ABS/PATH/TO/CICADA';
bidsRoot = '/ABS/PATH/TO/bids_data';
subject = '001';
session = '01';
runFolder = 'task-name';   % must match the completed CICADA run folder
icSelection = 'auto';          % 'manual' only for a completed Manual run

groupOutputDir = fullfile(bidsRoot, 'derivatives', 'cicada_group_qc', ...
    runFolder, 'CICADA_auto_ses_01'); % edit to the ACTUAL Group output dir
resultRoot = '/ABS/PATH/TO/NEW_INDEPENDENT_TASK_GLM_RESULTS';

% One-time code entry point. No recursive addpath/genpath needed.
addpath(fullfile(cicadaRepo, 'task_glm'), '-begin');

% 2. Resolve the exact producer-time recipe from this Group run.
% The recipe is per scan; do not select a different run's recipe.
scanDir = fullfile(bidsRoot, 'derivatives', 'cicada', ...
    ['sub-' subject], ['ses-' session], runFolder);
recipeName = sprintf('sub-%s_ses-%s_task-%s_matched_recipe.json', ...
    subject, session, runFolder);
recipeFile = fullfile(groupOutputDir, 'preparation_recipes', recipeName);

assert(isfolder(scanDir), 'Missing completed scan: %s', scanDir);
assert(isfile(recipeFile), ...
    ['Missing producer-attested Group recipe: %s. Do NOT substitute ', ...
     'a historical Group image or a guessed JSON.'], recipeFile);

% 3. Define contrast and covariance model BEFORE inspecting results.
% Condition names MUST exactly match the scan's hrf_conditions_raw.csv.
% These condition labels are foodpics examples, not hardcoded requirements.
cfg = struct();
cfg.scanDir = scanDir;
cfg.icSelection = icSelection;
cfg.inputMode = 'group_matched';
cfg.groupRecipe = recipeFile;
cfg.namedContrasts = struct('name', 'High_minus_Nonfood', ...
    'conditions', {{'highhedonic', 'nonfoods'}}, 'weights', [1, -1]);
cfg.driftModel = 'polynomial';
cfg.driftDegree = 2;
cfg.model = 'AR1';             % or 'AR2'; no automatic model selection
cfg.resultRoot = resultRoot;  % NEW analysis location, outside CICADA

% Optional: fit only scan-specific verified usable voxels, not entire
% broad ICA processing mask. The file must be in the SAME spatial grid.
% In older scans, check the canonical mask and retained voxel coverage.
cfg.analysisMaskFile = fullfile(scanDir, 'funcmask_constrained.nii.gz');

% 4. Safe first pass: verify lineage, timing, design, rank and mask.
% This is a true dry run: no BOLD fitting, no files written.
cfg.dryRun = true;
P = cicada_task_glm(cfg);
fprintf('Group task-GLM preflight: %s\n', P.status);

% 5. After preflight, deliberately opt in to a FULL first-level export.
% WARNING: this reads the full 4-D image; confirm sufficient RAM.
% It will refuse an existing output directory and will never produce
% calibrated significance maps. Set true ONLY when ready to fit/export.
doFullFitAndExport = false;
if doFullFitAndExport
    cfg.dryRun = false;
    cfg.runAnalysis = true;
    cfg.writeOutputs = true;
    cfg.allowGroupExport = true;
    cfg.maxVoxels = Inf;
    cfg.allowAllVoxels = true;
    cfg.allowLargeRun = true;
    R = cicada_task_glm(cfg);
    fprintf('Group first-level export: %s\n', R.status);
end
