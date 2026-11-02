function R = cicada_task_glm_reml_wholerun_readonly(S,opts)
% CICADA_TASK_GLM_REML_WHOLERUN_READONLY  Guarded scan-level REML/FGLS driver.
%
% Uses the unchanged, successfully tested numerical functions:
%   cicada_task_glm_reml_prepare.m
%   cicada_task_glm_reml_fit_batch.m
% and the preflight output S from the external FILM utility, called with
% cfg.execute=false. NO FSL command, image write, or other file write.
%
% IMPORTANT: Reads the pre-ICA-cleanup CICADA funcfile.nii.gz (NOT cleaned
% CICADA or Group output). May contain earlier despiking/preprocessing.
% This is fixed-ICA, model-conditional research code. Plug-in VARCOPE is
% provisional; do NOT convert to publication-grade p/z/t/CI inference.
%
% EXAMPLES (run after the existing FILM dry run):
%   W = cicada_task_glm_reml_wholerun_readonly(S);
%   W = cicada_task_glm_reml_wholerun_readonly(S, ...
%           struct('model','AR2','maxVoxels',256,'batchVoxels',64));
%   % EXACT prior 64-voxel sample (linear voxel indices):
%   W = cicada_task_glm_reml_wholerun_readonly(S135, ...
%           struct('model','AR2','voxelIndices',idx135,'batchVoxels',32));
%
% Deliberate defaults: AR1, up to 256 evenly spaced mask voxels, batches
% of 64, no full-brain invocation or output writing. To process >2048
% mask voxels (possibly expensive), explicitly pass allowLargeRun=true;
% maxVoxels=Inf additionally requires allowAllVoxels=true.
%
% MATLAB's niftiread loads the full compressed 4D image before batching.
% Thus the GLM is batchwise but input-file reading is NOT memory streaming.
% The output is IN MEMORY and spatially indexed by R.voxelIndices; it is
% not a NIfTI map until a separately validated export is implemented.

if nargin<2 || isempty(opts), opts=struct(); end
assert(isstruct(opts) && isscalar(opts),'opts must be a scalar struct');
opts = setDefault(opts,'model','AR1');
opts = setDefault(opts,'maxVoxels',256);
opts = setDefault(opts,'batchVoxels',64);
opts = setDefault(opts,'voxelIndices',[]);
opts = setDefault(opts,'allowLargeRun',false);
opts = setDefault(opts,'allowAllVoxels',false);
opts = setDefault(opts,'maxCacheBytes',256*1024^2);
opts = setDefault(opts,'returnBeta',false);
assert((ischar(opts.model) || (isstring(opts.model) && isscalar(opts.model))), ...
    'opts.model must be AR1 or AR2');
model=upper(char(string(opts.model)));
assert(ismember(model,{'AR1','AR2'}),'opts.model must be AR1 or AR2');
assert(isnumeric(opts.maxVoxels) && isscalar(opts.maxVoxels) && ...
    ~isnan(opts.maxVoxels) && opts.maxVoxels>=1 && ...
    (isinf(opts.maxVoxels) || opts.maxVoxels==floor(opts.maxVoxels)), ...
    'maxVoxels must be a positive integer or Inf');
assert(isnumeric(opts.batchVoxels) && isscalar(opts.batchVoxels) && ...
    isfinite(opts.batchVoxels) && opts.batchVoxels>=1 && ...
    opts.batchVoxels==floor(opts.batchVoxels) && opts.batchVoxels<=1024, ...
    'batchVoxels must be an integer between 1 and 1024');
assert(islogical(opts.allowLargeRun) && isscalar(opts.allowLargeRun), ...
    'allowLargeRun must be logical');
assert(islogical(opts.allowAllVoxels) && isscalar(opts.allowAllVoxels), ...
    'allowAllVoxels must be logical');
assert(islogical(opts.returnBeta) && isscalar(opts.returnBeta), ...
    'returnBeta must be logical');
assert(isnumeric(opts.maxCacheBytes) && isscalar(opts.maxCacheBytes) && ...
    isfinite(opts.maxCacheBytes) && opts.maxCacheBytes>0, ...
    'Invalid maxCacheBytes');
if isinf(opts.maxVoxels)
    assert(opts.allowAllVoxels, ...
        'maxVoxels=Inf requires opts.allowAllVoxels=true');
end

% Validate the ALREADY-DRY-RUN FILM preflight object. This function never
% calls that driver (and cannot accidentally trigger FILM execution).
assert(isstruct(S) && isscalar(S) && ...
    isfield(S,'status') && strcmp(S.status,'DRY_RUN') && ...
    isfield(S,'executed') && ~S.executed, ...
    'Supply S from a completed execute=false FILM preflight');
assert(isfield(S,'glm') && isstruct(S.glm) && ...
    all(isfield(S.glm,{'X0','C0','TR','sampleTimesSeconds'})) && ...
    isfield(S,'files') && all(isfield(S.files,{'original','mask'})), ...
    'Missing original-time model or source fields in S');
assert(isfield(S,'N') && isnumeric(S.N) && isscalar(S.N), ...
    'Missing volume count');
assert(isfile(S.files.original) && isfile(S.files.mask), ...
    'Original BOLD or scan mask no longer exists');
assert(size(S.glm.X0,1)==S.N && ...
    size(S.glm.C0,2)==size(S.glm.X0,2) && ...
    numel(S.glm.sampleTimesSeconds)==S.N, ...
    'Incorrect complete-model dimensions');

info=niftiinfo(S.files.original);
maskInfo=niftiinfo(S.files.mask);
assert(numel(info.ImageSize)>=4 && double(info.ImageSize(4))==S.N, ...
    'Original BOLD volume count disagrees with preflight');
spaceSize=double(info.ImageSize(1:3));
assert(numel(maskInfo.ImageSize)>=3 && ...
    isequal(double(maskInfo.ImageSize(1:3)),spaceSize), ...
    'Input mask spatial geometry mismatch');
assert(max(abs(double(info.Transform.T(:))- ...
    double(maskInfo.Transform.T(:))))<1e-3, ...
    'Input BOLD and mask affines differ');
assert(numel(info.PixelDimensions)>=4 && ...
    abs(double(info.PixelDimensions(4))-S.glm.TR)<1e-5, ...
    'Input BOLD TR changed since preflight');

mask=(niftiread(S.files.mask)>0);
assert(numel(mask)==prod(spaceSize), 'Mask image size mismatch');
maskIndex=find(mask(:));
assert(~isempty(maskIndex),'Empty CICADA mask');

if ~isempty(opts.voxelIndices)
    candidateIdx=double(opts.voxelIndices(:));
    assert(all(isfinite(candidateIdx)) && ...
        all(candidateIdx==floor(candidateIdx)) && ...
        all(candidateIdx>=1) && all(candidateIdx<=prod(spaceSize)), ...
        'voxelIndices must contain 1-based linear spatial voxel indices');
    assert(numel(unique(candidateIdx))==numel(candidateIdx), ...
        'Duplicate voxel indices are not supported');
    assert(all(mask(candidateIdx)), ...
        'At least one requested voxel is outside the CICADA mask');
    selection='explicit_voxel_indices';
else
    if isinf(opts.maxVoxels)
        candidateIdx=double(maskIndex);
    else
        nPick=min(numel(maskIndex),opts.maxVoxels);
        pick=round(linspace(1,numel(maskIndex),nPick));
        candidateIdx=double(maskIndex(pick));
    end
    selection='evenly_spaced_in_flattened_mask_order';
end
assert(~isempty(candidateIdx),'No selected voxels');
if numel(candidateIdx)>2048
    assert(opts.allowLargeRun, ...
        ['More than 2048 voxels requires allowLargeRun=true. ', ...
         'A full AR2 scan may be computationally expensive.']);
end
assert(numel(candidateIdx)<=prod(spaceSize),'Invalid voxel count');

% PREPARE ONCE: no covariance factorization is repeated per voxel batch.
cacheOpts=struct('model',model, ...
    'maxBatchVoxels',opts.batchVoxels, ...
    'maxCacheBytes',opts.maxCacheBytes);
t0=tic;
cache=cicada_task_glm_reml_prepare(S.glm.X0,S.glm.C0,cacheOpts);

% MATLAB niftiread decompresses/loads the full 4D BOLD; reshape is a
% logical view of the spatial voxel x time ordering (column-major).
B=niftiread(S.files.original);
assert(size(B,4)==S.N && ...
    (size(B,1)*size(B,2)*size(B,3))==prod(spaceSize), ...
    'Loaded BOLD dimensions differ from header');
Yall=reshape(B,[],S.N);
clear B

K=size(S.glm.C0,1);
V=numel(candidateIdx);
R=struct();
R.status='READ_ONLY_PROVISIONAL_FIXED_ICA';
R.model=model;
R.scanDir=S.scanDir;
R.sourceBOLD=S.files.original;
R.maskFile=S.files.mask;
R.noiseSelectionFile=S.files.noise;
R.rawHRFFile=S.files.rawHRF;
R.cleanedProvenanceJSON=S.files.cleanedJSON;
R.icSelection=S.icSelection;
R.driftModel=S.driftModel;
R.evNames=S.glm.evNames;
R.contrastMatrix=S.glm.C0;
R.X0Used = S.glm.X0;
R.batchVoxels=opts.batchVoxels;
R.maskVoxelCount=numel(maskIndex);
R.voxelIndices=candidateIdx(:);
R.selection=selection;
R.nRequestedVoxels=V;
R.valid=false(V,1);
R.nonfinite=false(V,1);
R.constant=false(V,1);
R.N=S.N;
R.TR=S.glm.TR;
R.imageSize=spaceSize;
R.imageAffine=double(info.Transform.T);
R.numRegressors=size(S.glm.X0,2);
R.residualDimension=cache.residualDimension;
R.nGridCandidates=cache.nGridCandidates;
R.approxCacheBytes=cache.approxCacheBytes;
R.contrastNames=S.contrasts;
R.taskTimingSource=S.taskTimingSource;
R.firstSampleSeconds=S.firstSampleSeconds;
R.r1Grid=cache.r1Grid;
R.k2Grid=cache.k2Grid;
R.cope=nan(K,V);
R.varcopePlugIn=nan(K,V);
R.r1Hat=nan(1,V);
R.kappa2Hat=nan(1,V);
R.remlProfileScore=nan(1,V);
R.residualScale2=nan(1,V);
R.r1AtGridBoundary=false(1,V);
R.kappa2AtGridBoundary=false(1,V);
if opts.returnBeta
    assert(V<=2048, ...
        'returnBeta=true limited to <=2048 voxels to bound memory');
    R.beta=nan(cache.P,V);
end

fprintf('CICADA original-BOLD REML | %s | %d candidates | batch=%d\n', ...
    model,V,opts.batchVoxels);
fprintf('Precomputed grid=%d, cache approximately %.1f MiB\n', ...
    cache.nGridCandidates,cache.approxCacheBytes/1024^2);
fprintf('READ ONLY: image data held in memory, no FSL/ICA/NIfTI output.\n');

for first=1:opts.batchVoxels:V
    last=min(V,first+opts.batchVoxels-1);
    pos=first:last;
    Y=double(Yall(candidateIdx(pos),:)).';
    finite=all(isfinite(Y),1);
    variable=(max(Y,[],1)>min(Y,[],1));
    valid=finite & variable;
    R.nonfinite(pos)=~finite(:);
    R.constant(pos)=(finite & ~variable)';
    R.valid(pos)=valid(:);
    if ~any(valid), continue; end
    cols=pos(valid);
    fit=cicada_task_glm_reml_fit_batch(Y(:,valid),cache);
    R.cope(:,cols)=fit.cope;
    R.varcopePlugIn(:,cols)=fit.varcopePlugIn;
    R.r1Hat(cols)=fit.r1Hat;
    R.kappa2Hat(cols)=fit.kappa2Hat;
    R.remlProfileScore(cols)=fit.remlProfileScore;
    R.residualScale2(cols)=fit.residualScale2;
    R.r1AtGridBoundary(cols)=fit.r1AtGridBoundary;
    R.kappa2AtGridBoundary(cols)=fit.kappa2AtGridBoundary;
    if opts.returnBeta
        R.beta(:,cols)=fit.beta;
    end
end
clear Yall
R.nValidVoxels=nnz(R.valid);
R.nInvalidNonfinite=nnz(R.nonfinite);
R.nInvalidConstant=nnz(R.constant);
R.elapsedSeconds=toc(t0);
R.notes=['Results are indexed by voxelIndices (1-based linear locations ' ...
    'in original spatial volume). Invalid voxels have NaN results. ' ...
    'All COPE/VARCOPE values are exploratory fixed-ICA conditional GLS; ' ...
    'VARCOPE excludes estimation of covariance parameters and adaptive ' ...
    'ICA/classification uncertainty. No valid p-value calibration ' ...
    'is claimed. No files were written.'];
fprintf('%s | valid=%d/%d | nonfinite=%d | constant=%d | %.2f s\n', ...
    model,R.nValidVoxels,V,R.nInvalidNonfinite, ...
    R.nInvalidConstant,R.elapsedSeconds);
end

function opts=setDefault(opts,name,value)
if ~isfield(opts,name) || isempty(opts.(name)), opts.(name)=value; end
end
