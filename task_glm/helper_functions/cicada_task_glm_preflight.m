function S = cicada_task_glm_preflight(cfg)
% CICADA_TASK_GLM_PREFLIGHT  Read-only original-BOLD task design preflight.
% Builds the existing task HRF + drift + selected ICA-noise design, named
% contrasts, run metadata and full-rank original-time model for REML/FGLS.
% Does NOT run FILM, fit BOLD, create directories or write any files.
%
% This routine preserves the task-design mathematical recipe and output
% fields of the November 13, 2026 legacy FILM dry-run preflight. Its default
% execution mode is always nonexecuting. cfg.execute=true is rejected.
% The former experimental FILM pilot is no longer included in the
% task_glm companion; this preflight never invokes a FILM comparator.
%
% For task GLM, use the public cicada_task_glm.m entry point.

if nargin < 1 || isempty(cfg), cfg = struct(); end
assert(isstruct(cfg) && isscalar(cfg), 'cfg must be a scalar struct');

cfg = default(cfg, 'bidsRoot', '');
cfg = default(cfg, 'subject', '');
cfg = default(cfg, 'session', '');
cfg = default(cfg, 'runFolder', '');
cfg = default(cfg, 'scanDir', '');
cfg = default(cfg, 'expectedTR', []);
cfg = default(cfg, 'contrastNames', {});
cfg = default(cfg, 'contrastTaskWeights', []);
cfg = default(cfg, 'namedContrasts', []);
cfg = default(cfg, 'outputRoot', '');
cfg = default(cfg, 'driftDegree', 2);
cfg = default(cfg, 'driftModel', 'polynomial');
cfg = default(cfg, 'highPassHz', 0.008);
cfg = default(cfg, 'execute', false);
cfg = default(cfg, 'outputDir', '');
cfg = default(cfg, 'referenceDesignMat', '');
cfg = default(cfg, 'icSelection', 'auto');
cfg = default(cfg, 'analysisMaskFile', '');
assert(isscalar(cfg.execute) && ...
    (islogical(cfg.execute) || ...
    (isnumeric(cfg.execute) && ...
    isfinite(cfg.execute) && ...
    ismember(cfg.execute,[0 1]))), ...
    'execute must be true/false or 0/1');
assert(~logical(cfg.execute), ...
    'Design preflight never executes FILM or a task GLM; set execute=false');
assert(ismember(cfg.driftDegree, [0 1 2]), 'Use drift degree 0/1/2');
assert(ischar(cfg.subject) || isstring(cfg.subject));
assert(ischar(cfg.session) || isstring(cfg.session));
assert(ischar(cfg.runFolder) || isstring(cfg.runFolder));
selection = lower(char(string(cfg.icSelection)));
assert(ismember(selection, {'auto','manual'}), ...
    'icSelection must be auto or manual');

driftModel = lower(char(string(cfg.driftModel)));

assert(ismember(driftModel, {'polynomial','cosine'}), ...
    'driftModel must be polynomial or cosine');

% Allow the user to specify any completed CICADA scan directory.
% Otherwise preserve the existing BIDS-based default.

if ~isempty(cfg.scanDir)
    assert((ischar(cfg.scanDir) && isrow(cfg.scanDir)) || ...
           (isstring(cfg.scanDir) && isscalar(cfg.scanDir)), ...
        'cfg.scanDir must be a path string');

    scanDir = char(cfg.scanDir);
else
    assert(~isempty(cfg.bidsRoot) && ...
           ~isempty(cfg.subject) && ...
           ~isempty(cfg.session) && ...
           ~isempty(cfg.runFolder), ...
        ['Provide cfg.scanDir, or specify bidsRoot, ', ...
         'subject, session, and runFolder']);

    scanDir = fullfile(cfg.bidsRoot, 'derivatives', 'cicada', ...
        char(cfg.subject), char(cfg.session), char(cfg.runFolder));
end

assert(isfolder(scanDir), ...
    'CICADA scan directory not found: %s', scanDir);
files.original = fullfile(scanDir,'funcfile.nii.gz');

% Discover the completed Auto CICADA cleaned image and JSON.
% Require exactly one match to avoid silently choosing a wrong run.

cleanedDir = fullfile(scanDir, 'cleaned');

J = dir(fullfile(cleanedDir, ...
    sprintf('*CICADA_%s_nonagg_bold.json', selection)));

% Ignore macOS AppleDouble metadata files.
J = J(~startsWith(string({J.name}), "._"));

assert(numel(J) == 1, ...
    'Expected one %s CICADA cleaned-BOLD JSON in %s; found %d: %s', ...
    selection, cleanedDir, numel(J), ...
    strjoin({J.name}, ', '));

files.cleanedJSON = fullfile(J(1).folder, J(1).name);
files.cleaned = regexprep(files.cleanedJSON, ...
    '\.json$', '.nii.gz');

assert(isfile(files.cleaned), ...
    'Missing cleaned BOLD corresponding to %s', files.cleanedJSON);
files.mix = fullfile(scanDir,'melodic','melodic_mix');
files.noise = fullfile(scanDir, ...
    sprintf('ic_%s_selection',selection), ...
    sprintf('%s_noise_dist_ICs.csv',selection));
files.rawHRF = fullfile(scanDir,'hrf_conditions_raw.csv');
% The Group/ICA processing mask is scientifically fixed by the producer
% recipe. Never replace it when selecting voxels for the downstream GLM.
files.mask = fullfile(scanDir,'funcmask.nii.gz');
files.analysisMask = files.mask;
analysisMaskSelection = 'processing_default';
if ~isempty(cfg.analysisMaskFile)
    assert((ischar(cfg.analysisMaskFile) && isrow(cfg.analysisMaskFile)) || ...
        (isstring(cfg.analysisMaskFile) && isscalar(cfg.analysisMaskFile)), ...
        'cfg.analysisMaskFile must be a file path');
    requested = char(cfg.analysisMaskFile);
    assert(~isempty(strtrim(requested)) && ...
        java.io.File(requested).isAbsolute(), ...
        'cfg.analysisMaskFile must be an absolute path');
    files.analysisMask = requested;
    analysisMaskSelection = 'explicit_binary_subset';
end
fields = {'original','mix','noise','rawHRF','mask', ...
    'analysisMask','cleanedJSON'};
for k = 1:numel(fields)
    assert(isfile(files.(fields{k})), 'Missing %s: %s',fields{k},files.(fields{k}));
end

% Only inspect the original BOLD header and scan's existing mask geometry.
info = niftiinfo(files.original);
maskInfo = niftiinfo(files.mask);
assert(numel(info.ImageSize) >= 4, 'Input BOLD is not 4D');
N = double(info.ImageSize(4));
assert(isequal(double(info.ImageSize(1:3)),double(maskInfo.ImageSize(1:3))), ...
    'BOLD/mask grids differ');
assert(max(abs(double(info.Transform.T(:))-double(maskInfo.Transform.T(:)))) < 1e-3, ...
    'BOLD/mask affines differ');

% Validate the independent GLM sampling mask before any fitting or writing.
% No atlas resampling, threshold guessing, or modifying the ICA mask.
analysisInfo = niftiinfo(files.analysisMask);
assert(numel(analysisInfo.ImageSize)==3 && ...
    isequal(double(analysisInfo.ImageSize(:)'), ...
    double(info.ImageSize(1:3))), ...
    'Analysis mask must be a 3D image on the exact BOLD grid');
assert(max(abs(double(info.Transform.T(:))- ...
    double(analysisInfo.Transform.T(:)))) < 1e-3, ...
    'Analysis mask affine differs from BOLD');
processingMaskData = double(niftiread(files.mask));
analysisMaskData = double(niftiread(files.analysisMask));
assert(all(isfinite(analysisMaskData(:))) && ...
    all(analysisMaskData(:)==0 | analysisMaskData(:)==1), ...
    'Analysis mask must contain only finite binary values 0 or 1');
assert(any(analysisMaskData(:)>0), 'Analysis mask is empty');
assert(~any(analysisMaskData(:)>0 & processingMaskData(:)<=0), ...
    'Analysis mask extends outside CICADA processing mask');
processingMaskVoxelCount = nnz(processingMaskData>0);
analysisMaskVoxelCount = nnz(analysisMaskData>0);
clear processingMaskData analysisMaskData

if ~isempty(cfg.expectedTR)
    assert(numel(info.PixelDimensions)>=4, 'Missing NIfTI TR');
    assert(abs(double(info.PixelDimensions(4))-cfg.expectedTR)<1e-5, ...
        'Unexpected TR: inspect the input and run metadata');
end

% Verify the CICADA HRF sampling reference used for this scan.
meta = jsondecode(fileread(files.cleanedJSON));

assert(isfield(meta,'TaskHRF'), ...
    'Missing TaskHRF provenance in cleaned-BOLD JSON');

timing = meta.TaskHRF;

assert(isfield(timing,'timing_source') && ...
       isfield(timing,'sample_times_seconds'), ...
    'Incomplete CICADA task timing provenance');

sampleTimes = double(timing.sample_times_seconds(:));
TR = double(info.PixelDimensions(4));

if strcmp(driftModel,'cosine')
    assert(isnumeric(cfg.highPassHz) && ...
        isscalar(cfg.highPassHz) && ...
        isfinite(cfg.highPassHz) && ...
        cfg.highPassHz > 0 && ...
        cfg.highPassHz < 1/(2*TR), ...
        'Invalid cosine high-pass frequency');
end

assert(numel(sampleTimes)==N && ...
       all(isfinite(sampleTimes)), ...
    'Invalid task sampling times');

assert(all(abs(diff(sampleTimes)-TR)<1e-5), ...
    'Task sampling times inconsistent with BOLD TR');

M = readmatrix(files.mix,'FileType','text');
noise = readmatrix(files.noise);
noise = double(noise(:)');
T = readtable(files.rawHRF,'VariableNamingRule','preserve');
H = double(T{:,:});
labels = T.Properties.VariableNames;
assert(size(M,1)==N && size(H,1)==N, 'Time-point mismatch');
assert(all(isfinite(M(:))) && all(isfinite(H(:))), 'Non-finite matrix');
assert(all(isfinite(noise)) && all(noise==floor(noise)) && ...
    all(noise>=1) && all(noise<=size(M,2)) && ...
    numel(unique(noise))==numel(noise), 'Invalid/duplicate ICA noise IDs');

assert(isfield(meta,'ICSelection') && ...
       isfield(meta.ICSelection,'NoiseICs'), ...
    'Missing ICA selection provenance');

sidecarNoise = double(meta.ICSelection.NoiseICs(:)');

assert(isequal(sort(sidecarNoise),sort(noise)), ...
    'Selected noise ICs disagree with cleaned-BOLD provenance');

nTask = size(H,2);
assert(nTask >= 1, 'No task conditions in raw HRF table');

assert(numel(unique(labels)) == nTask, ...
    'Duplicate raw-HRF condition names');


assert(rank(M)==size(M,2), 'MELODIC matrix rank deficient');

% Center the task, drift and selected-ICA columns while retaining an
% explicit intercept in the REML design returned in S.glm.X0.

X = H;
evNames = labels;

switch driftModel

    case 'polynomial'
        t = linspace(-1,1,N)';

        if cfg.driftDegree >= 1
            X = [X t];
            evNames{end+1} = 'linear_drift';
        end

        if cfg.driftDegree >= 2
            X = [X t.^2-mean(t.^2)];
            evNames{end+1} = 'quadratic_drift';
        end

    case 'cosine'
        % DCT-II drift basis, matching Nilearn's convention.
        nDCT = floor(2*N*TR*cfg.highPassHz);
        sampleIndex = (0:N-1)';

        for k = 1:nDCT
            d = sqrt(2/N) * cos( ...
                pi*(sampleIndex+0.5)*k/N);

            X = [X d]; %#ok<AGROW>
            evNames{end+1} = sprintf( ...
                'cosine_drift_%02d',k); %#ok<AGROW>
        end
end

nDrift = size(X,2)-nTask;

X = [X M(:,noise)];
evNames = [evNames, arrayfun(@(i)sprintf('noise_IC_%d',i),noise, ...
    'UniformOutput',false)];
X = X - mean(X,1);
p = size(X,2);
r = rank(X);
assert(r==p && rank([ones(N,1) X])==p+1, ...
    'Joint design not full rank; model requires review');
assert(N-p-1>0, 'No residual dimensions');


if isempty(cfg.namedContrasts)
    assert(~isempty(cfg.contrastTaskWeights) && ...
           ~isempty(cfg.contrastNames), ...
        'Specify namedContrasts or explicit contrast weights and names');

    % Original numeric-contrast interface remains available.
    CW = double(cfg.contrastTaskWeights);
    contrastNames = cfg.contrastNames;
else
    % Resolve requested contrasts by task-condition name.
    defs = cfg.namedContrasts;
    assert(isstruct(defs), ...
        'namedContrasts must be a struct array');

    CW = zeros(numel(defs), nTask);
    contrastNames = cell(1, numel(defs));

    for j = 1:numel(defs)
        d = defs(j);

        assert(all(isfield(d, ...
            {'name','conditions','weights'})), ...
            'Each contrast needs name, conditions, and weights');

        conds = cellstr(string(d.conditions));
        w = double(d.weights(:)');

        [found, loc] = ismember(conds, labels);

        assert(all(found), ...
            'Contrast %d refers to an unknown condition', j);
        assert(numel(unique(loc)) == numel(loc), ...
            'Contrast %d contains duplicate conditions', j);
        assert(numel(w) == numel(conds) && ...
            all(isfinite(w)) && any(w ~= 0), ...
            'Invalid weights for contrast %d', j);

        CW(j,loc) = w;
        contrastNames{j} = char(string(d.name));

        assert(~isempty(strtrim(contrastNames{j})), ...
            'Empty contrast name');
    end

    assert(numel(unique(contrastNames)) == ...
        numel(contrastNames), 'Duplicate contrast names');
end

assert(size(CW,2) == nTask && ...
    size(CW,1) == numel(contrastNames), ...
    'Contrast definitions do not match the HRF table');

assert(all(isfinite(CW(:))) && ...
    all(any(CW ~= 0,2)), 'Invalid contrast');

C = [CW zeros(size(CW,1), p-nTask)];

% Preflight design diagnostics.
% Compare each task contrast with versus without ICA noise EVs.
% These are unweighted design-only diagnostics, not REML standard errors.

nBase = nTask + nDrift;

[~, Rfull] = qr(X, 0);
[~, Rbase] = qr(X(:,1:nBase), 0);

fullFactor = sum((Rfull' \ C').^2, 1)';
baseFactor = sum((Rbase' \ C(:,1:nBase)').^2, 1)';

contrastVIF = fullFactor ./ baseFactor;

% Condition number after normalizing EV column magnitudes.
colNorm = sqrt(sum(X.^2,1));
assert(all(colNorm > 0), 'Zero-valued design column');

scaledX = X ./ colNorm;
singularValues = svd(scaledX, 'econ');
scaledCondition = singularValues(1) / singularValues(end);

assert(all(isfinite(contrastVIF)) && ...
       all(contrastVIF >= 1-1e-8), ...
       'Unexpected contrast variance-inflation diagnostic');

S = struct();
S.scanDir = scanDir;
S.files = files;
S.N = N;
S.numICs = size(M,2);
S.numNoiseICs = numel(noise);
S.numWaves = p;
S.designRank = r;
S.residualDimension = N-p-1;
S.evs = evNames;
S.contrasts = contrastNames;
S.taskContrastWeights = CW;
S.contrastVIF = contrastVIF;
S.scaledDesignCondition = scaledCondition;
S.icSelection = selection;
S.driftModel = driftModel;
S.nDriftRegressors = nDrift;
S.contrastDesignVarianceFactor = fullFactor;
% Complete original-time joint task and nuisance design for REML/FGLS.

S.glm = struct();
S.glm.X0 = [ones(N,1), X];
S.glm.C0 = [zeros(size(C,1),1), C];

S.glm.evNames = [{'intercept'}, evNames];
S.glm.taskNames = labels;
S.glm.noiseICs = noise;

S.glm.TR = TR;
S.glm.sampleTimesSeconds = sampleTimes;
S.glm.rank = r + 1;
S.glm.residualDimension = N - r - 1;

% Resolve analysis output directory without mislabeling custom scans.
defaultOutputDir = fullfile(cfg.outputRoot, ...
    char(cfg.subject), char(cfg.session), ...
    [char(cfg.runFolder) '_' selection '_joint_reml']);

if ~isempty(cfg.outputDir)

    % Explicit output location takes precedence.
    assert((ischar(cfg.outputDir) && isrow(cfg.outputDir)) || ...
           (isstring(cfg.outputDir) && isscalar(cfg.outputDir)), ...
        'cfg.outputDir must be a path string');

    S.outputDir = char(cfg.outputDir);

else

    if ~isempty(cfg.scanDir)
        % Without an explicit output directory, the labels supplied
        % in cfg must match the selected scan directory.
        [parentDir, runLabel] = fileparts(scanDir);
        [parentDir, sessionLabel] = fileparts(parentDir);
        [~, subjectLabel] = fileparts(parentDir);

        assert(strcmp(runLabel, char(cfg.runFolder)) && ...
               strcmp(sessionLabel, char(cfg.session)) && ...
               strcmp(subjectLabel, char(cfg.subject)), ...
            ['Custom scanDir does not match the subject/session/run ', ...
             'labels. Set cfg.outputDir explicitly to avoid ', ...
             'mislabeling results.']);
    end
    
    assert(~isempty(cfg.outputRoot) && ...
       ~isempty(cfg.subject) && ...
       ~isempty(cfg.session) && ...
       ~isempty(cfg.runFolder), ...
    ['Specify cfg.outputDir, or provide outputRoot ', ...
     'and subject/session/runFolder labels']);
    S.outputDir = defaultOutputDir;
end

% Prevent writes into CICADA-managed derivative trees.
outReal = char(java.io.File(S.outputDir).getCanonicalPath());
scanReal = char(java.io.File(scanDir).getCanonicalPath());

managedRoot = fileparts(fileparts(fileparts(scanReal)));
[~, managedName] = fileparts(managedRoot);

if strcmp(managedName, 'cicada')
    protectedRoots = {managedRoot, ...
        fullfile(fileparts(managedRoot), 'cicada_group_qc')};
else
    protectedRoots = {scanReal};
end

for k = 1:numel(protectedRoots)
    rootReal = char(java.io.File( ...
        protectedRoots{k}).getCanonicalPath());

    assert(~strcmp(outReal, rootReal) && ...
           ~startsWith(outReal, [rootReal filesep]), ...
        'Output directory is inside CICADA-managed data');
end

S.executed = false;
S.status = 'DRY_RUN';
S.referenceDesignMaxAbsError = NaN;
S.maskVoxelCount = analysisMaskVoxelCount;
S.processingMaskVoxelCount = processingMaskVoxelCount;
S.analysisMaskSelection = analysisMaskSelection;
S.taskTimingSource = char(string(timing.timing_source));
S.firstSampleSeconds = sampleTimes(1);

% Optional reproducibility comparison with a saved VEST design matrix.
if ~isempty(cfg.referenceDesignMat)
    assert(isfile(cfg.referenceDesignMat),'Reference design missing');
    Xref = readVestMatrix(cfg.referenceDesignMat);
    assert(isequal(size(Xref),size(X)), 'Reference design size mismatch');
    S.referenceDesignMaxAbsError = max(abs(Xref(:)-X(:)));
    fprintf('Max difference vs reference FILM design: %.9g\n', ...
        S.referenceDesignMaxAbsError);
    assert(S.referenceDesignMaxAbsError < 1e-8, ...
        'Generated design differs from the reference design');
end

fprintf('CICADA task GLM | design preflight only | no fitting or writes\n');
fprintf('%s | N=%d, EVs=%d, rank=%d, residual=%d\n', ...
    scanDir,N,p,r,S.residualDimension);
fprintf('GLM analysis mask voxels: %d (processing mask: %d)\n', ...
    S.maskVoxelCount,S.processingMaskVoxelCount);
fprintf('Output location (not created): %s\n',S.outputDir);
fprintf('PREFLIGHT ONLY: original inputs unchanged, no FILM invoked.\n');

end

function cfg = default(cfg,name,value)
if ~isfield(cfg,name) || isempty(cfg.(name))
    cfg.(name)=value;
end
end

function X=readVestMatrix(path)
fid=fopen(path,'r');
assert(fid>0,'Cannot read reference design');
cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
lines={};
while ~feof(fid)
    x=fgetl(fid);
    if ischar(x),lines{end+1}=x; end %#ok<AGROW>
end
m=find(strcmp(strtrim(lines),'/Matrix'),1);
assert(~isempty(m),'Invalid VEST design: no /Matrix');
rowLines=lines(m+1:end);
rowLines=rowLines(~cellfun(@isempty,rowLines));
rows=cellfun(@(s)sscanf(s,'%f')',rowLines,'UniformOutput',false);
assert(~isempty(rows) && all(cellfun(@numel,rows)==numel(rows{1})), ...
    'Invalid VEST rows');
X=vertcat(rows{:});
end
