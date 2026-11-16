function R = cicada_roi_timeseries(cfg)
% CICADA_ROI_TIMESERIES  Public, atlas-independent BOLD ROI extraction.
%
% R = cicada_roi_timeseries(cfg)
%
% cfg.boldFile: 4D BOLD NIfTI (prefer unsmoothed for anatomical ROIs)
% cfg.roiFiles: one or multiple 3D/4D ROI NIfTIs
% cfg.roiTypes: corresponding 'labels', 'binary', or 'weighted'
% cfg.qualityMaskFile: aligned analysis-quality mask; alternatively,
%   cfg.scanDir resolves <scanDir>/funcmask_constrained.nii.gz
% cfg.dryRun=true: geometry and ROI coverage preflight, no time-series read
% cfg.writeOutputs=true: guarded CSV/JSON export into a NEW output path
%
% Uses the independently tested v0.2 extraction implementation. For
% authenticated CSV/JSON export, it temporarily resolves CICADA's existing
% SHA-256 helper. It never runs task GLM, ICA, or FSL and does not leave
% helper directories on the user's MATLAB path.
%
% Status: experimental source integration; NOT a release validation.

assert(nargin==1 && isstruct(cfg) && isscalar(cfg), ...
    'Provide one scalar cfg struct');
roiDir=fileparts(mfilename('fullpath'));
helperDir=fullfile(roiDir,'helper_functions');
impl=fullfile(helperDir,'cicada_roi_timeseries_v02.m');
assert(isfile(impl), ...
    'Missing ROI implementation: %s',impl);

repoRoot=fileparts(roiDir);
hashDir=fullfile(repoRoot,'task_glm','helper_functions');
hashFile=fullfile(hashDir,'cicada_task_glm_sha256.m');
needsHash=isfield(cfg,'writeOutputs') && isequal(cfg.writeOutputs,true);
if needsHash
    assert(isfile(hashFile), ...
        'Missing CICADA SHA-256 helper required for ROI export: %s',hashFile);
end

priorPath=path;
restorePath=onCleanup(@()path(priorPath)); %#ok<NASGU>
if isfile(hashFile)
    addpath(hashDir,'-begin');
end
addpath(helperDir,'-begin');
assert(strcmp(which('cicada_roi_timeseries_v02'),impl), ...
    'Incorrect ROI implementation resolved on MATLAB path');
if needsHash
    assert(strcmp(which('cicada_task_glm_sha256'),hashFile), ...
        'Incorrect CICADA SHA-256 helper resolved on MATLAB path');
end
R=cicada_roi_timeseries_v02(cfg);
end
