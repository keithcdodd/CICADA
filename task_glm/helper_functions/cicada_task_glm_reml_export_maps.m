function P = cicada_task_glm_reml_export_maps(W,S,opts)
% CICADA_TASK_GLM_REML_EXPORT_MAPS  Guarded exploratory map export.
%
% P = cicada_task_glm_reml_export_maps(W,S,opts)
%   W: read-only output of cicada_task_glm_reml_wholerun_readonly
%   S: SAME scan's execute=false FILM preflight (provides full X0/C0)
%
% No calculation, FILM, denoising, or source changes. No writes by default.
% Explicit writeOutputs=true creates ONLY a new output directory, never
% overwrites, and stages files until their spatial/data checks pass.
%
% Default writeOutputs=false: preview and validate only.
% Sparse W requires allowSparseTrial=true AND uses SPARSE_TRIAL labels.
% A sparse map is NOT a whole-brain result; NaN marks untested voxels.
% Outputs: sampled/valid masks; per-contrast COPE and provisional VARCOPE;
% r1 map (plus kappa2 for AR2); analysis_model.mat; provenance.json.
% No t/z/p maps; NO calibrated inference is claimed.
%
% MATLAB niftiwrite is used with the ORIGINAL 3D mask header so affine,
% voxel geometry, and orientation can be verified after writing.
%
% Example preview (ZERO writes):
%  o=struct('outputDir','/Users/me/task_glm_trial_sparse64');
%  P=cicada_task_glm_reml_export_maps(W135,S135,o);
% Explicit one-time sparse export:
%  o.writeOutputs=true; o.allowSparseTrial=true;
%  P=cicada_task_glm_reml_export_maps(W135,S135,o);
%
% An independently verified original-time mean/covariance model remains
% scientifically PROVISIONAL because ICA estimation/labels were adaptive.

if nargin<3 || isempty(opts), opts=struct(); end
assert(isstruct(opts) && isscalar(opts),'opts must be a scalar struct');
opts=def(opts,'outputDir','');
opts=def(opts,'writeOutputs',false);
opts=def(opts,'allowSparseTrial',false);
assert(isscalar(opts.writeOutputs) && islogical(opts.writeOutputs), ...
    'writeOutputs must be true or false (logical)');
assert(isscalar(opts.allowSparseTrial) && islogical(opts.allowSparseTrial), ...
    'allowSparseTrial must be true or false (logical)');
assert(ischar(opts.outputDir) || ...
    (isstring(opts.outputDir) && isscalar(opts.outputDir)), ...
    'outputDir must be a path');
assert(~isempty(strtrim(char(opts.outputDir))), ...
    'Specify an independent outputDir');
assert(isstruct(W) && isscalar(W) && ...
    isfield(W,'status') && ...
    strcmp(W.status,'READ_ONLY_PROVISIONAL_FIXED_ICA'), ...
    'W must be output from the read-only whole-run driver');
assert(isstruct(S) && isscalar(S) && isfield(S,'status') && ...
    strcmp(S.status,'DRY_RUN') && ~S.executed, ...
    'S must be the matching execute=false FILM preflight');
assert(all(isfield(W,{'sourceBOLD','maskFile','scanDir','imageSize', ...
    'imageAffine','N','TR','model','voxelIndices','valid', ...
    'maskVoxelCount','nRequestedVoxels','nValidVoxels', ...
    'cope','varcopePlugIn','r1Hat','kappa2Hat', ...
    'r1AtGridBoundary','kappa2AtGridBoundary', ...
    'residualScale2','contrastMatrix','contrastNames', ...
    'r1Grid','k2Grid','nGridCandidates'})), ...
    'W is missing required fields');
assert(isfield(S,'glm') && ...
    all(isfield(S.glm,{'X0','C0','TR','sampleTimesSeconds'})), ...
    'S is missing the original-time design');
assert(isfield(S,'files') && ...
    all(isfield(S.files,{'original','mask','mix','noise', ...
    'rawHRF','cleanedJSON'})), ...
    'S is missing required source paths');

% Verify matching lineage and all available design identities.
assert(strcmp(canon(S.scanDir),canon(W.scanDir)), ...
    'W and S use different scan directories');
assert(strcmp(canon(S.files.original),canon(W.sourceBOLD)), ...
    'W and S refer to different original BOLD files');
assert(strcmp(canon(S.files.mask),canon(W.maskFile)), ...
    'W and S refer to different masks');
assert(W.N==S.N && W.TR==S.glm.TR, ...
    'W and S temporal dimensions differ');
assert(W.numRegressors==size(S.glm.X0,2) && ...
    W.residualDimension==S.glm.residualDimension, ...
    'W and S model dimensions differ');
assert(isequal(W.contrastMatrix,S.glm.C0), ...
    'W and S contrast matrices differ');
assert(strcmpi(W.model,'AR1') || strcmpi(W.model,'AR2'), ...
    'Unrecognized covariance model');
assert(isfile(W.maskFile) && isfile(W.sourceBOLD), ...
    'Referenced BOLD or mask no longer exists');

maskInfo=niftiinfo(W.maskFile);
boldInfo=niftiinfo(W.sourceBOLD);
spaceSize=double(W.imageSize(:)');
assert(numel(spaceSize)==3 && all(spaceSize>0) && ...
    all(spaceSize==floor(spaceSize)), ...
    'Invalid spatial dimensions');
assert(isequal(spaceSize,double(maskInfo.ImageSize(1:3))) && ...
    isequal(spaceSize,double(boldInfo.ImageSize(1:3))), ...
    'Spatial dimensions changed since W was computed');
maskAffine=double(maskInfo.Transform.T);
assert(isequal(size(maskAffine),[4 4]) && ...
    max(abs(maskAffine(:)-W.imageAffine(:)))<1e-3 && ...
    max(abs(maskAffine(:)-double(boldInfo.Transform.T(:))))<1e-3, ...
    'Spatial affine changed since W was computed');
Nvox=prod(spaceSize);
idx=double(W.voxelIndices(:));
V=numel(idx);
valid=logical(W.valid(:));
assert(V>0 && W.nRequestedVoxels==V && ...
    nnz(valid)==W.nValidVoxels && ...
    all(isfinite(idx)) && all(idx>=1) && ...
    all(idx<=Nvox) && all(idx==floor(idx)) && ...
    numel(unique(idx))==V, ...
    'Invalid spatial indexing or voxel counts');
assert(all(isfield(W,{'numRegressors','residualDimension'})), ...
    'W lacks design dimension fields');
mask=niftiread(W.maskFile)>0;
assert(numel(mask)==Nvox && ...
    nnz(mask)==W.maskVoxelCount && all(mask(idx)), ...
    'Scan mask disagrees with W');
K=size(W.contrastMatrix,1);
assert(K>=1 && numel(W.contrastNames)==K && ...
    isequal(size(W.cope),[K V]) && ...
    isequal(size(W.varcopePlugIn),[K V]), ...
    'Contrast map array dimensions do not match W');
assert(isfield(W,'X0Used') && ...
    isequaln(W.X0Used,S.glm.X0), ...
    'The original-time design changed between fit and export');
assert(all(isfinite(W.cope(:,valid)),'all') && ...
    all(isfinite(W.varcopePlugIn(:,valid)) & ...
    W.varcopePlugIn(:,valid)>0,'all'), ...
    'Valid voxels must have finite COPE and positive VARCOPE');
assert(all(isfinite(W.r1Hat(valid))) && ...
    all(isfinite(W.kappa2Hat(valid))) && ...
    all(isfinite(W.residualScale2(valid))) && ...
    all(W.residualScale2(valid)>0), ...
    'Invalid covariance estimates or residual scales');
assert(all(isnan(W.cope(:,~valid)),'all') && ...
    all(isnan(W.varcopePlugIn(:,~valid)),'all'), ...
    'Invalid fit positions must not contain contrast estimates');

% Sparse means not every processing-mask voxel was attempted.
% An attempted but invalid voxel is not an untested voxel.
isSparse = (V ~= W.maskVoxelCount);
if isSparse && opts.writeOutputs
    assert(opts.allowSparseTrial, ...
        ['Sparse or incomplete W: writeOutputs=true requires ', ...
        'allowSparseTrial=true. Sparse outputs are QC-only.']);
end
scope='COMPLETE_MASK_ATTEMPT';
if isSparse, scope='SPARSE_TRIAL'; end

% Canonical-path ancestry guards, evaluated even during dry-run.
outDir=canon(char(opts.outputDir));
scanDir=canon(S.scanDir);
managedRoot=canon(fileparts(fileparts(fileparts(scanDir))));
[~,rootLabel]=fileparts(managedRoot);
protected={scanDir,canon(W.maskFile),canon(W.sourceBOLD)};
if strcmp(rootLabel,'cicada')
    protected={scanDir,managedRoot, ...
        canon(fullfile(fileparts(managedRoot),'cicada_group_qc'))};
end
for q=1:numel(protected)
    r=protected{q};
    assert(~strcmp(outDir,r) && ...
        ~startsWith(outDir,[r filesep]), ...
        'Output directory is inside a protected CICADA source tree');
end
assert(~isfolder(outDir) && ~isfile(outDir), ...
    'Output already exists (no overwrite): %s',outDir);

% Before a WRITE, compare the current scientific inputs and source code
% against hashes captured before the fit. Previews remain fast/no-write.
% Older in-memory v0.2 fits cannot provide this integrity guarantee.
if opts.writeOutputs
    assert(isfield(W,'fileFingerprints') && ...
        isstruct(W.fileFingerprints), ...
        ['This fit has no pre-fit SHA256 fingerprints. ', ...
         'Use the matching historical exporter for historical results; ', ...
         'new exports require a v0.3 companion fit.']);
    currentFingerprints=cicada_task_glm_fingerprints(S);
    assert(isequaln(W.fileFingerprints,currentFingerprints), ...
        'Input or companion source content changed since fit; refusing export');
end
parent = fileparts(outDir);
assert(~isempty(parent), ...
    'An absolute output directory is required');

P=struct();
P.status='PREVIEW_ONLY';
P.outputDir=outDir;
P.scope=scope;
P.model=upper(W.model);
P.sparseTrial=isSparse;
P.numContrasts=K;
P.nRequested=V;
P.nValid=nnz(valid);
P.maskVoxelCount=W.maskVoxelCount;
filePrefix='';
if isSparse, filePrefix='sparseTrial_'; end
P.wouldWrite={'analysis_model.mat','provenance.json', ...
    [filePrefix 'desc-sampled_mask.nii.gz'], ...
    [filePrefix 'desc-valid_mask.nii.gz'], ...
    [filePrefix 'desc-r1Hat.nii.gz']};
if strcmpi(W.model,'AR2')
    P.wouldWrite{end+1}=[filePrefix 'desc-kappa2Hat.nii.gz'];
end
for k=1:K
    P.wouldWrite{end+1}=sprintf('%scontrast-%02d_cope.nii.gz',filePrefix,k);
    P.wouldWrite{end+1}=sprintf('%scontrast-%02d_varcopePlugIn.nii.gz',filePrefix,k);
end
fprintf('REML export preview | %s | %s | %d/%d valid, mask=%d\n', ...
    P.model,scope,P.nValid,V,W.maskVoxelCount);
fprintf('Destination (new directory only): %s\n',outDir);
fprintf('No inferential p/t/z maps; VARCOPE is provisional.\n');
if ~opts.writeOutputs
    fprintf('PREVIEW ONLY. No files written.\n');
    return
end
% Create the result parent only when writing is authorized.
% All output-path and input-integrity checks have already passed.
if ~isfolder(parent)
    [ok,msg] = mkdir(parent);
    assert(ok && isfolder(parent), ...
        'Could not create result directory %s: %s',parent,msg);
end

% Write into a unique sibling staging folder. Only promote after all
% files and voxel-grid/affine/data roundtrip checks have succeeded.
stageDir=[outDir '.partial_' char(java.util.UUID.randomUUID)];
assert(~isfolder(stageDir) && ~isfile(stageDir), ...
    'Unexpected staging-directory collision');
[ok,msg]=mkdir(stageDir);
assert(ok,'Cannot create staging directory: %s',msg);
try
    % The maskInfo carries the 3D spatial transform. Update datatype and
    % bit depth for floating-point statistical images, retaining geometry.
    floatInfo=maskInfo;
    floatInfo.Datatype='single';
    floatInfo.BitsPerPixel=32;
    for rf={'raw','Raw'}
        n=rf{1};
        if isfield(floatInfo,n)
            if isfield(floatInfo.(n),'datatype')
                floatInfo.(n).datatype=16;
            end
            if isfield(floatInfo.(n),'bitpix')
                floatInfo.(n).bitpix=32;
            end
        end
    end
    floatInfo.Description='CICADA provisional fixed-ICA REML/FGLS';
    maskWriteInfo=maskInfo;
    maskWriteInfo.Datatype='uint8';
    maskWriteInfo.BitsPerPixel=8;
    for rf={'raw','Raw'}
        n=rf{1};
        if isfield(maskWriteInfo,n)
            if isfield(maskWriteInfo.(n),'datatype')
                maskWriteInfo.(n).datatype=2;
            end
            if isfield(maskWriteInfo.(n),'bitpix')
                maskWriteInfo.(n).bitpix=8;
            end
        end
    end
    maskWriteInfo.Description='CICADA REML export coverage mask';

    sampledMask=zeros(Nvox,1,'uint8');
    validMask=zeros(Nvox,1,'uint8');
    sampledMask(idx)=1;
    validMask(idx(valid))=1;
    emit(reshape(sampledMask,spaceSize), ...
        [filePrefix 'desc-sampled_mask'],maskWriteInfo,stageDir,spaceSize,maskAffine);
    emit(reshape(validMask,spaceSize), ...
        [filePrefix 'desc-valid_mask'],maskWriteInfo,stageDir,spaceSize,maskAffine);
    emitFloat(W.r1Hat,[filePrefix 'desc-r1Hat']);
    if strcmpi(W.model,'AR2')
        emitFloat(W.kappa2Hat,[filePrefix 'desc-kappa2Hat']);
    end
    for k=1:K
        emitFloat(W.cope(k,:),sprintf('%scontrast-%02d_cope',filePrefix,k));
        emitFloat(W.varcopePlugIn(k,:), ...
            sprintf('%scontrast-%02d_varcopePlugIn',filePrefix,k));
    end

    % Save the EXACT numerical design used by the checked preflight S.
    X0=S.glm.X0; C0=S.glm.C0;
    evNames=S.glm.evNames; contrastNames=S.contrasts;
    sampleTimesSeconds=S.glm.sampleTimesSeconds;
    noiseICs=S.glm.noiseICs;
    r1Grid=W.r1Grid; k2Grid=W.k2Grid;
    save(fullfile(stageDir,'analysis_model.mat'),'X0','C0', ...
        'evNames','contrastNames','sampleTimesSeconds','noiseICs', ...
        'r1Grid','k2Grid','-v7');

    manifest=struct();
    manifest.schema='cicada_external_reml_spatial_export_v0.2';
    manifest.status='EXPLORATORY_PROVISIONAL_FIXED_ICA';
    manifest.scope=scope;
    manifest.model=P.model;
    manifest.sourceBOLD=W.sourceBOLD;
    manifest.sourceMask=W.maskFile;
    manifest.scanDir=S.scanDir;
    manifest.selection=S.icSelection;
    manifest.inputFiles=S.files;
    manifest.inputFileMetadata=inputStamps(S.files);
    manifest.companionVersion=W.fileFingerprints.companionVersion;
    manifest.matlabVersion=W.fileFingerprints.matlabVersion;
    manifest.contentFingerprints=W.fileFingerprints;
    manifest.volumeCount=W.N;
    manifest.TR=W.TR;
    manifest.imageSize=spaceSize;
    manifest.voxelAffine=maskAffine;
    manifest.nRegressors=W.numRegressors;
    manifest.nominalResidualDimension=W.residualDimension;
    manifest.nGridCandidates=W.nGridCandidates;
    manifest.r1Grid=W.r1Grid;
    manifest.k2Grid=W.k2Grid;
    manifest.contrastNames=S.contrasts;
    manifest.contrastWeights=S.glm.C0;
    manifest.nMaskVoxels=W.maskVoxelCount;
    manifest.nAttemptedVoxels=V;
    manifest.nValidVoxels=nnz(valid);
    manifest.nNonfiniteVoxels=nnz(W.nonfinite);
    manifest.nConstantVoxels=nnz(W.constant);
    manifest.nR1GridBoundaryVoxels=nnz(W.r1AtGridBoundary);
    manifest.nKappa2GridBoundaryVoxels=nnz(W.kappa2AtGridBoundary);
    manifest.taskTimingSource=S.taskTimingSource;
    manifest.firstSampleSeconds=S.firstSampleSeconds;
    manifest.generatedAtLocal=char(datetime('now','Format','yyyy-MM-dd''T''HH:mm:ss'));
    manifest.dataConventions=['All floating-point maps are NaN where ', ...
        'unfitted, invalid, or untested. desc-sampled_mask identifies ', ...
        'attempted voxels; desc-valid_mask identifies successful fits. ', ...
        'Sparse pilot maps must never be treated as whole-brain results.'];
    manifest.statisticalLimitations=['Original pre-ICA-cleanup BOLD with ', ...
        'selected ICA noise ICs jointly modeled; fixed-ICA conditional ', ...
        'stationary AR covariance. VARCOPE is a plug-in estimator and ', ...
        'not validated for real-data p/t/z inference, adaptive ICA, ', ...
        'model misspecification or spatial multiple testing.'];
    manifest.provenanceLimitations=['SHA256 fingerprints of the original ', ...
        'BOLD, mask, MELODIC matrix, selected noise ICs, raw HRF and ', ...
        'cleaned JSON were computed before model fitting and checked ', ...
        'before export. Numerical design and covariance grids are ', ...
        'stored in analysis_model.mat. No cryptographic guarantee ', ...
        'protects against a file being changed during the fit itself. ', ...
        'The cleaned BOLD image is not part of the numerical input.'];
    manifest.generatedFiles=P.wouldWrite;
    fid=fopen(fullfile(stageDir,'provenance.json'),'w');
    assert(fid>0,'Unable to write provenance.json');
    fprintf(fid,'%s\n',jsonencode(manifest,'PrettyPrint',true));
    fclose(fid);

    assert(~isfolder(outDir) && ~isfile(outDir), ...
        'Output directory appeared during staging; refusing overwrite');
    [ok,msg]=movefile(stageDir,outDir);
    assert(ok,'Failed to finalize staged export: %s',msg);
catch ME
    if isfolder(stageDir)
        rmdir(stageDir,'s'); % Only our newly created staging directory.
    end
    rethrow(ME)
end
P.status='EXPORTED_PROVISIONAL_NO_OVERWRITE';
P.exportedFiles=P.wouldWrite;
fprintf('EXPORTED and reread-verified %d files to %s\n', ...
    numel(P.exportedFiles),outDir);

    function emitFloat(values,stem)
        values=double(values(:));
        assert(numel(values)==V,'Map data count mismatch');
        assert(all(isfinite(single(values(valid)))) && ...
            all(isfinite(values(valid))), ...
            'Map values exceed Float32 finite precision');
        map=nan(Nvox,1,'single');
        map(idx(valid))=single(values(valid));
        emit(reshape(map,spaceSize),stem,floatInfo, ...
            stageDir,spaceSize,maskAffine);
    end
end

function emit(vol,stem,info,folder,spaceSize,referenceAffine)
% Filename is an extension-free STEM: Compressed=true generates .nii.gz.
filePath=fullfile(folder,[stem '.nii.gz']);
assert(~isfile(filePath),'Unexpected map overwrite');
niftiwrite(vol,fullfile(folder,stem),info,'Compressed',true);
assert(isfile(filePath),'niftiwrite did not produce expected .nii.gz');
outInfo=niftiinfo(filePath);
assert(isequal(double(outInfo.ImageSize(1:3)),spaceSize), ...
    'Exported map spatial shape changed');
assert(max(abs(double(outInfo.Transform.T(:))-referenceAffine(:)))<1e-3, ...
    'Exported map affine differs from source mask');
back=niftiread(filePath);
assert(isequaln(back,vol), ...
    'Exported map value roundtrip mismatch');
end

function names=inputStamps(files)
keys=fieldnames(files);
names=struct();
for q=1:numel(keys)
    f=char(files.(keys{q}));
    if ~isfile(f),continue;end
    d=dir(f);
    names.(keys{q})=struct('path',f,'bytes',d.bytes, ...
        'modifiedLocal',d.date);
end
end

function s=canon(path)
s=char(java.io.File(char(path)).getCanonicalPath());
end

function opts=def(opts,name,value)
if ~isfield(opts,name) || isempty(opts.(name))
    opts.(name)=value;
end
end
