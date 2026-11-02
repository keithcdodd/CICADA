function R = cicada_task_glm(cfg)
% CICADA_TASK_GLM  Optional, standalone CICADA task-contrast companion.
%
% R = cicada_task_glm(cfg)
%
% Single-run, original-BOLD joint task + selected ICA-noise GLM.
% Reads ALREADY COMPLETED Auto/Manual CICADA scan products; never calls
% CICADA denoising or changes its managed derivatives. No Group CICADA.
% Requires CICADA hrf_conditions_raw.csv and matching TaskHRF provenance.
% No inference-calibrated t/z/p/CI output is produced.
%
% Normal use: one call validates, fits the full processing mask,
% and exports to a NEW, separate output directory (never overwrites).
% Supply cfg.resultRoot for automatic subject/session/run/mode naming,
% or cfg.outputDir for an explicitly chosen result directory.
% Optional cfg.dryRun=true: validate/design only, with no writes.
% Legacy cfg.runAnalysis and cfg.writeOutputs options continue working.
%
% Required: cfg.scanDir (or bidsRoot+subject+session+runFolder),
%           cfg.resultRoot or cfg.outputDir,
%           cfg.namedContrasts (or contrastNames+contrastTaskWeights).
% Model: cfg.model='AR1' (default) or 'AR2'.
% Input: original CICADA funcfile.nii.gz; NOT cleaned BOLD.
% Full-mask processing loads the entire 4D image into memory before
% statistical batching; check available memory before use.
%
% This is a guarded experimental companion, not a released inferential
% analysis or an automatic CICADA processing stage.

if nargin<1 || isempty(cfg), cfg=struct(); end
assert(isstruct(cfg) && isscalar(cfg), ...
    'cfg must be a scalar configuration struct');

% Friendly default: run and save on a normal call. Retain the old
% development-mode switches when explicitly supplied by an older script.
legacyPhase=isfield(cfg,'runAnalysis') || isfield(cfg,'writeOutputs');
cfg=opt(cfg,'dryRun',false);
assert(islogical(cfg.dryRun) && isscalar(cfg.dryRun), ...
    'cfg.dryRun must be logical true or false');
if cfg.dryRun
    assert(~(isfield(cfg,'runAnalysis') && isequal(cfg.runAnalysis,true)) && ...
           ~(isfield(cfg,'writeOutputs') && isequal(cfg.writeOutputs,true)), ...
           'dryRun=true conflicts with explicit runAnalysis/writeOutputs=true');
    cfg.runAnalysis=false;
    cfg.writeOutputs=false;
elseif ~legacyPhase
    cfg.runAnalysis=true;
    cfg.writeOutputs=true;
    cfg=opt(cfg,'maxVoxels',Inf);
    cfg=opt(cfg,'allowAllVoxels',true);
    cfg=opt(cfg,'allowLargeRun',true);
end

cfg=opt(cfg,'runAnalysis',false);
cfg=opt(cfg,'writeOutputs',false);
cfg=opt(cfg,'model','AR1');
cfg=opt(cfg,'resultRoot','');
cfg=opt(cfg,'icSelection','auto');
cfg=opt(cfg,'driftModel','cosine');
cfg=opt(cfg,'highPassHz',0.008);
cfg=opt(cfg,'maxVoxels',256);
cfg=opt(cfg,'voxelIndices',[]);
cfg=opt(cfg,'batchVoxels',64);
cfg=opt(cfg,'allowAllVoxels',false);
cfg=opt(cfg,'allowLargeRun',false);
cfg=opt(cfg,'allowSparseTrial',false);
cfg=opt(cfg,'returnBeta',false);
cfg=opt(cfg,'maxCacheBytes',256*1024^2);

for f={'runAnalysis','writeOutputs','allowAllVoxels','allowLargeRun', ...
        'allowSparseTrial','returnBeta'}
    name=f{1};
    assert(islogical(cfg.(name)) && isscalar(cfg.(name)), ...
        'cfg.%s must be scalar logical true/false',name);
end
assert(~cfg.writeOutputs || cfg.runAnalysis, ...
    'writeOutputs=true requires runAnalysis=true');
assert(ischar(cfg.model) || ...
    (isstring(cfg.model) && isscalar(cfg.model)), ...
    'cfg.model must be AR1 or AR2');
cfg.model=upper(char(string(cfg.model)));
assert(ismember(cfg.model,{'AR1','AR2'}), ...
    'cfg.model must be AR1 or AR2');

% Do not accept legacy FILM execution through this wrapper. The proven
% FILM comparator remains separately available for explicitly requested
% use; it has different execution/output semantics.
if isfield(cfg,'execute')
    assert(isequal(cfg.execute,false) || isequal(cfg.execute,0), ...
        ['Do not set cfg.execute=true here; use ', ...
         'cfg.runAnalysis and cfg.writeOutputs']);
end

% Compute scan-specific outputDir from an optional resultRoot. The root
% need not exist: the exporter creates it only for an authorized write.
if ~isempty(cfg.resultRoot)
    assert(ischar(cfg.resultRoot) || ...
        (isstring(cfg.resultRoot) && isscalar(cfg.resultRoot)), ...
        'cfg.resultRoot must be a path string');
    root=char(cfg.resultRoot);
    assert(~isempty(strtrim(root)) && java.io.File(root).isAbsolute(), ...
        'cfg.resultRoot must be an absolute directory path');
    assert((~isfield(cfg,'outputDir') || isempty(cfg.outputDir)) && ...
           (~isfield(cfg,'outputRoot') || isempty(cfg.outputRoot)), ...
        'Specify resultRoot OR outputDir/outputRoot, not both');
    if isfield(cfg,'scanDir') && ~isempty(cfg.scanDir)
        [sparent,runLabel]=fileparts(char(cfg.scanDir));
        [sparent,sesLabel]=fileparts(sparent);
        [~,subLabel]=fileparts(sparent);
    else
        needed={'subject','session','runFolder'};
        for ii=1:numel(needed)
            assert(isfield(cfg,needed{ii}) && ~isempty(cfg.(needed{ii})), ...
                'resultRoot requires subject/session/runFolder or scanDir');
        end
        subLabel=char(cfg.subject);
        sesLabel=char(cfg.session);
        runLabel=char(cfg.runFolder);
    end
    selectionLabel=lower(char(string(cfg.icSelection)));
    folderLabels={subLabel,sesLabel,runLabel,selectionLabel};
    for ii=1:numel(folderLabels)
        assert(~isempty(regexp(folderLabels{ii}, ...
            '^[A-Za-z0-9][A-Za-z0-9_-]*$','once')), ...
            'Unsafe result-folder label: %s',folderLabels{ii});
    end
    outputLabel=strjoin(folderLabels,'_');
    cfg.outputDir=fullfile(root,outputLabel);
    % Do not accidentally put analyses inside the CICADA source tree.
    repoRoot=char(java.io.File( ...
        fileparts(fileparts(mfilename('fullpath')))).getCanonicalPath());
    rootReal=char(java.io.File(root).getCanonicalPath());
    assert(~strcmp(rootReal,repoRoot) && ...
           ~startsWith(rootReal,[repoRoot filesep]), ...
        'resultRoot must be outside the CICADA source repository');
end
assert(isfield(cfg,'outputDir') && ...
    (ischar(cfg.outputDir) || ...
     (isstring(cfg.outputDir) && isscalar(cfg.outputDir))) && ...
    ~isempty(strtrim(char(cfg.outputDir))), ...
    'Set cfg.resultRoot or cfg.outputDir outside CICADA derivatives');

% Prevent misspelled execution options from silently doing a dry run.
known={'runAnalysis','writeOutputs','dryRun','resultRoot', ...
    'model','icSelection', ...
    'driftModel','highPassHz','maxVoxels','voxelIndices', ...
    'batchVoxels','allowAllVoxels','allowLargeRun','allowSparseTrial', ...
    'returnBeta','maxCacheBytes','execute','outputDir','scanDir', ...
    'bidsRoot','subject','session','runFolder','expectedTR', ...
    'contrastNames','contrastTaskWeights','namedContrasts', ...
    'outputRoot','driftDegree','referenceDesignMat', ...
    'filmThreshold','filmAcfSpatialSmoothing','filmAcfMaskSize'};
unknown=setdiff(fieldnames(cfg),known);
assert(isempty(unknown), ...
    'Unsupported cfg field(s): %s',strjoin(unknown,', '));

% Keep implementation functions hidden from the normal user workflow.
% Add this companion's own helpers only for the duration of this call;
% restore the caller's original MATLAB path even if an error occurs.
helperDir = fullfile(fileparts(mfilename('fullpath')), 'helper_functions');
assert(isfolder(helperDir), ...
    'Missing task_glm/helper_functions. Keep the companion folder intact');
originalPath = path;
restoreHelperPath = onCleanup(@()path(originalPath)); %#ok<NASGU>
addpath(helperDir, '-begin');

% Reuse the validated preflight verbatim; force no execution.
preCfg=cfg;
preCfg.execute=false;
if cfg.runAnalysis
    % Save verbose preflight messages; do not print the helper's historical
    % FILM dry-run messages during an ordinary full analysis.
    preflightLog=evalc('S=cicada_task_glm_option1_film_pilot(preCfg);');
    fprintf('CICADA task GLM | %s | %d volumes | %s\n', ...
        S.scanDir,S.N,cfg.model);
else
    preflightLog='';
    S=cicada_task_glm_option1_film_pilot(preCfg);
end
assert(strcmp(S.status,'DRY_RUN') && ~S.executed, ...
    'Preflight unexpectedly attempted execution');

R=struct();
R.status='PREFLIGHT_ONLY';
R.method='ORIGINAL_BOLD_JOINT_REML_FGLS';
R.model=cfg.model;
R.preflight=S;
R.preflightLog=preflightLog;
R.fit=[];
R.export=[];
R.notes=['No t/z/p maps or calibrated confidence intervals. ', ...
    'ICA classification is data-adaptive and covariance estimation may ', ...
    'produce anti-conservative plug-in uncertainty.'];

if ~cfg.runAnalysis
    fprintf('OPTIONAL TASK GLM: PREFLIGHT ONLY. No fit or files written.\n');
    return
end

% Fail early before expensive image loading/fitting when output already
% exists, or a sparse export was requested without explicit permission.
assert(~isfolder(S.outputDir) && ~isfile(S.outputDir), ...
    'Output already exists; refusing overwrite: %s',S.outputDir);
if cfg.writeOutputs && ~cfg.allowSparseTrial
    % Explicitly selected voxels are always treated as a sparse trial
    % unless they actually enumerate the whole mask.
    explicit=~isempty(cfg.voxelIndices);
    if explicit
        fullRequested=numel(cfg.voxelIndices)==S.maskVoxelCount;
    else
        fullRequested=isinf(cfg.maxVoxels) || ...
            cfg.maxVoxels>=S.maskVoxelCount;
    end
    assert(fullRequested, ...
        ['Writing sampled maps requires cfg.allowSparseTrial=true. ', ...
         'For a complete mask set maxVoxels=Inf and both allow flags.']);
end

fitOpts=struct('model',cfg.model, ...
    'maxVoxels',cfg.maxVoxels, ...
    'voxelIndices',cfg.voxelIndices, ...
    'batchVoxels',cfg.batchVoxels, ...
    'allowAllVoxels',cfg.allowAllVoxels, ...
    'allowLargeRun',cfg.allowLargeRun, ...
    'maxCacheBytes',cfg.maxCacheBytes, ...
    'returnBeta',cfg.returnBeta);

% Hash all numerically used scan inputs and companion source before fitting.
% This records actual file content, not merely paths or timestamps.
fileFingerprints=cicada_task_glm_fingerprints(S);
W=cicada_task_glm_reml_wholerun_readonly(S,fitOpts);
W.fileFingerprints=fileFingerprints;
R.fit=W;
R.status='IN_MEMORY_FIT_PROVISIONAL';

% Run the same guarded exporter first in PREVIEW mode. This checks
% source lineage, exact design, geometry, coverage, and path safety.
exportOpts=struct('outputDir',S.outputDir, ...
    'writeOutputs',false, ...
    'allowSparseTrial',cfg.allowSparseTrial);
P=cicada_task_glm_reml_export_maps(W,S,exportOpts);
R.export=P;

if cfg.writeOutputs
    exportOpts.writeOutputs=true;
    P=cicada_task_glm_reml_export_maps(W,S,exportOpts);
    R.export=P;
    R.status=P.status;
end

fprintf('OPTIONAL TASK GLM: %s\n',R.status);
end

function cfg=opt(cfg,name,value)
if ~isfield(cfg,name) || isempty(cfg.(name))
    cfg.(name)=value;
end
end
