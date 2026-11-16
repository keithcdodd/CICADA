function R = cicada_roi_timeseries_v02(cfg)
% CICADA_ROI_TIMESERIES_V02  Atlas-independent (optionally weighted) BOLD ROI means.
% External experimental helper; NOT yet integrated into CICADA source.
%
% Required:
%   cfg.boldFile: 4D NIfTI (time in fourth dimension).
%   cfg.roiFiles: one ROI/atlas NIfTI path or cell array of paths.
% Optional:
%   cfg.roiTypes: one type per file, one of 'labels', 'binary', 'weighted'.
%     Default is 'labels' for backward compatibility with v01.
%     'labels': one 3D nonnegative integer parcellation, each positive
%         integer becomes a separate ROI. A binary map also works here.
%     'binary': one 3D 0/1 mask or 4D stack of independently named 0/1 ROIs.
%     'weighted': a 3D map or 4D stack of finite nonnegative weights.
%         Each 3D map / stack volume becomes ONE ROI, even if empty.
%         Overlap between weighted maps in 4D is allowed.
%         Values need not sum to 1 across maps or be <= 1.
%         Explicit 'weighted' prevents fractional maps from being
%         mistaken for label atlases.
%   cfg.qualityMaskFile: binary 3D mask; OR cfg.scanDir to resolve
%     <scanDir>/funcmask_constrained.nii.gz. No silent fallback.
%   cfg.allowNoQualityMask=true explicitly disables quality masking.
%   cfg.minCoverageFraction (0..1, default 0). For weighted ROIs it
%     applies to retained weight mass; otherwise retained voxel count.
%   cfg.minValidVoxels (positive integer, default 1).
%   cfg.dryRun, cfg.writeOutputs, cfg.outputDir (NEW directory only).
%
% Strict native-grid checks; NO implicit resampling or registration.
% Each voxel must be finite, nonzero and temporally nonconstant for
% extraction. Weights are renormalized after quality and signal screens.
% Reports both voxel and weight-mass coverage, weighted effective N.
% Time series with inadequate coverage are NaN rather than zero.
% No inferential statistics, FSL calls or writes to existing derivatives.

assert(nargin==1 && isstruct(cfg) && isscalar(cfg), ...
    'Provide one scalar cfg struct');
cfg=local_def(cfg,'scanDir','');
cfg=local_def(cfg,'roiTypes',{});
cfg=local_def(cfg,'qualityMaskFile','');
cfg=local_def(cfg,'allowNoQualityMask',false);
cfg=local_def(cfg,'minCoverageFraction',0);
cfg=local_def(cfg,'minValidVoxels',1);
cfg=local_def(cfg,'dryRun',false);
cfg=local_def(cfg,'writeOutputs',false);
cfg=local_def(cfg,'outputDir','');
assert(isfield(cfg,'boldFile') && isfield(cfg,'roiFiles'), ...
    'cfg.boldFile and cfg.roiFiles are required');
assert(local_flag(cfg.dryRun) && local_flag(cfg.writeOutputs) && ...
       local_flag(cfg.allowNoQualityMask), ...
    'dryRun, writeOutputs and allowNoQualityMask must be logical flags');
assert(~(cfg.dryRun && cfg.writeOutputs), ...
    'A dry-run cannot write outputs');
assert(isnumeric(cfg.minCoverageFraction) && ...
    isscalar(cfg.minCoverageFraction) && ...
    isfinite(cfg.minCoverageFraction) && ...
    cfg.minCoverageFraction>=0 && cfg.minCoverageFraction<=1, ...
    'minCoverageFraction must be in [0,1]');
assert(isnumeric(cfg.minValidVoxels) && isscalar(cfg.minValidVoxels) && ...
    isfinite(cfg.minValidVoxels) && cfg.minValidVoxels>=1 && ...
    cfg.minValidVoxels==floor(cfg.minValidVoxels), ...
    'minValidVoxels must be a positive integer');

boldFile=local_path(cfg.boldFile);
assert(isfile(boldFile),'Missing BOLD: %s',boldFile);
roiFiles=cfg.roiFiles;
if ischar(roiFiles) || (isstring(roiFiles) && isscalar(roiFiles))
    roiFiles={char(roiFiles)};
elseif isstring(roiFiles)
    roiFiles=cellstr(roiFiles(:));
end
assert(iscell(roiFiles) && ~isempty(roiFiles), ...
    'roiFiles must be a nonempty filename or cell array of filenames');
roiFiles=cellfun(@local_path,roiFiles,'UniformOutput',false);
assert(all(cellfun(@isfile,roiFiles)), ...
    'At least one requested ROI NIfTI is missing');
roiTypes=cfg.roiTypes;
if isempty(roiTypes)
    roiTypes=repmat({'labels'},size(roiFiles));
elseif ischar(roiTypes) || (isstring(roiTypes) && isscalar(roiTypes))
    roiTypes={char(roiTypes)};
elseif isstring(roiTypes)
    roiTypes=cellstr(roiTypes(:));
end
assert(iscell(roiTypes) && numel(roiTypes)==numel(roiFiles), ...
    'roiTypes must specify one type per roiFiles entry');
for fi=1:numel(roiTypes)
    type=roiTypes{fi};
    assert((ischar(type) && isrow(type)) || ...
        (isstring(type) && isscalar(type)), ...
        'Invalid roiTypes element');
    type=lower(strtrim(char(type)));
    assert(ismember(type,{'labels','binary','weighted'}), ...
        'roiTypes must be labels, binary, or weighted');
    roiTypes{fi}=type;
end

qualityFile='';
if ~isempty(cfg.qualityMaskFile)
    qualityFile=local_path(cfg.qualityMaskFile);
elseif ~cfg.allowNoQualityMask
    assert(~isempty(cfg.scanDir), ...
        ['Provide qualityMaskFile, or scanDir to use ', ...
         'funcmask_constrained.nii.gz, or explicitly set ', ...
         'allowNoQualityMask=true']);
    qualityFile=fullfile(local_path(cfg.scanDir),'funcmask_constrained.nii.gz');
end
if ~isempty(qualityFile)
    assert(isfile(qualityFile),'Missing quality mask: %s',qualityFile);
end

boldInfo=niftiinfo(boldFile);
assert(numel(boldInfo.ImageSize)==4 && ...
    numel(boldInfo.PixelDimensions)>=4 && ...
    all(double(boldInfo.ImageSize)>0), 'BOLD must be 4D with a TR');
dims=double(boldInfo.ImageSize(1:3));
N=double(boldInfo.ImageSize(4));
TR=double(boldInfo.PixelDimensions(4));
assert(isfinite(TR) && TR>0,'Invalid BOLD TR');
nVox=prod(dims);

if isempty(qualityFile)
    quality=true(nVox,1);
    qualitySelection='explicitly_disabled';
else
    qInfo=niftiinfo(qualityFile);
    local_same_grid(boldInfo,qInfo,qualityFile,false);
    q=double(niftiread(qualityFile));
    assert(all(isfinite(q(:))) && all(q(:)==0 | q(:)==1), ...
        'Quality mask must contain only 0/1 finite values');
    quality=reshape(q>0,[],1);
    qualitySelection='verified_binary_mask';
end

labels={};
roiFileColumn={};
labelIDs=[];
volumeIndices=[];
roiKinds={};
roiIndices={};
roiWeights={};
for fi=1:numel(roiFiles)
    info=niftiinfo(roiFiles{fi});
    kind=roiTypes{fi};
    if strcmp(kind,'labels') && numel(info.ImageSize)==4
        error('4D ROI stacks require explicit binary or weighted roiTypes');
    end
    local_same_grid(boldInfo,info,roiFiles{fi},~strcmp(kind,'labels'));
    sz=double(info.ImageSize);
    if numel(sz)==3
        nMaps=1;
    else
        nMaps=sz(4);
    end
    assert(nMaps>=1 && nMaps<=4096,'Invalid ROI stack length');
    values=double(niftiread(roiFiles{fi}));
    assert(numel(values)==nVox*nMaps, ...
        'ROI image data/header size mismatch: %s',roiFiles{fi});
    values=reshape(values,nVox,nMaps);
    assert(all(isfinite(values(:))) && all(values(:)>=0), ...
        'ROI values must be finite and nonnegative: %s',roiFiles{fi});
    if strcmp(kind,'labels')
        assert(nMaps==1 && all(values(:)==floor(values(:))), ...
            'Label atlas must be a 3D nonnegative integer map: %s',roiFiles{fi});
        present=unique(values(values>0));
        assert(~isempty(present),'ROI file has no positive labels: %s',roiFiles{fi});
        assert(numel(present)<=4096,'Too many ROI labels in one atlas');
        for li=1:numel(present)
            code=present(li);
            ids=find(values(:,1)==code);
            [labels,roiFileColumn,labelIDs,volumeIndices,roiKinds,roiIndices,roiWeights] = ...
                local_append_roi(labels,roiFileColumn,labelIDs,volumeIndices, ...
                roiKinds,roiIndices,roiWeights,roiFiles{fi},ids, ...
                ones(numel(ids),1),code,1,'labels');
        end
    else
        if strcmp(kind,'binary')
            assert(all(values(:)==0 | values(:)==1), ...
                'Binary ROI must contain only 0 and 1: %s',roiFiles{fi});
        end
        for vi=1:nMaps
            ids=find(values(:,vi)>0);
            if strcmp(kind,'binary')
                w=ones(numel(ids),1);
                code=1;
            else
                w=values(ids,vi);
                code=NaN;
            end
            [labels,roiFileColumn,labelIDs,volumeIndices,roiKinds,roiIndices,roiWeights] = ...
                local_append_roi(labels,roiFileColumn,labelIDs,volumeIndices, ...
                roiKinds,roiIndices,roiWeights,roiFiles{fi},ids, ...
                w,code,vi,kind);
        end
    end
end
K=numel(labels);
assert(K<=4096,'Too many total ROIs');
atlasN=zeros(K,1);
qualityN=zeros(K,1);
validN=NaN(K,1);
atlasMass=zeros(K,1);
qualityMass=zeros(K,1);
validMass=NaN(K,1);
effectiveN=NaN(K,1);
status=strings(K,1);
for ri=1:K
    ids=roiIndices{ri};
    w=roiWeights{ri};
    atlasN(ri)=numel(ids);
    qualityN(ri)=nnz(quality(ids));
    atlasMass(ri)=sum(w);
    qualityMass(ri)=sum(w(quality(ids)));
    assert(isfinite(atlasMass(ri)) && isfinite(qualityMass(ri)), ...
        'Weight-mass sum overflow in ROI %d',ri);
    if atlasN(ri)==0
        status(ri)="empty_roi";
    else
        status(ri)="preflight_only";
    end
end
qualityFrac=local_safe_fraction(qualityN,atlasN);
validFrac=NaN(K,1);
qualityMassFrac=local_safe_fraction(qualityMass,atlasMass);
validMassFrac=NaN(K,1);
coverageMetric=strings(K,1);
for ri=1:K
    if strcmp(roiKinds{ri},'weighted')
        coverageMetric(ri)="weight_mass";
    else
        coverageMetric(ri)="voxel_fraction";
    end
end
series=NaN(N,K);
if ~cfg.dryRun
    Y=niftiread(boldFile);
    assert(numel(Y)==nVox*N,'BOLD data/header size mismatch');
    Y=reshape(Y,nVox,N);
    for ri=1:K
        ids=roiIndices{ri};
        if isempty(ids)
            validN(ri)=0;
            validMass(ri)=0;
            effectiveN(ri)=0;
            status(ri)="empty_roi";
            continue
        end
        w=roiWeights{ri};
        keep=quality(ids);
        ids=ids(keep);
        w=w(keep);
        if isempty(ids)
            validN(ri)=0;
            validMass(ri)=0;
            effectiveN(ri)=0;
            status(ri)="no_quality_coverage";
            continue
        end
        subset=double(Y(ids,:));
        validRow=all(isfinite(subset),2) & ...
            all(subset~=0,2) & ...
            any(subset~=subset(:,1),2);
        validN(ri)=nnz(validRow);
        validW=w(validRow);
        validMass(ri)=sum(validW);
        if validMass(ri)>0
            effectiveN(ri)=validMass(ri)^2/sum(validW.^2);
        else
            effectiveN(ri)=0;
        end
        if validN(ri)==0
            status(ri)="no_valid_signal";
        else
            if strcmp(roiKinds{ri},'weighted')
                retainedFraction=validMass(ri)/atlasMass(ri);
            else
                retainedFraction=validN(ri)/atlasN(ri);
            end
            if validN(ri)<cfg.minValidVoxels || ...
                    retainedFraction<cfg.minCoverageFraction
                status(ri)="low_coverage";
            else
                % Positive weights renormalized over accepted voxels.
                normalizedW=validW/validMass(ri);
                series(:,ri)=subset(validRow,:).'*normalizedW;
                status(ri)="ok";
            end
        end
    end
    validFrac=local_safe_fraction(validN,atlasN);
    validMassFrac=local_safe_fraction(validMass,atlasMass);
end
coverage=table(string(labels),string(roiFileColumn),labelIDs, ...
    atlasN,qualityN,validN,qualityFrac,validFrac,status, ...
    string(roiKinds),volumeIndices,atlasMass,qualityMass,validMass, ...
    qualityMassFrac,validMassFrac,effectiveN,coverageMetric, ...
    'VariableNames',{'roi','roi_file','label_id','atlas_voxels', ...
    'quality_voxels','valid_voxels','quality_fraction', ...
    'valid_fraction','status','roi_type','volume_index', ...
    'atlas_weight_mass','quality_weight_mass','valid_weight_mass', ...
    'quality_weight_fraction','valid_weight_fraction', ...
    'effective_valid_voxels','coverage_metric'});

R=struct();
R.status='PREFLIGHT_ONLY';
if ~cfg.dryRun, R.status='EXTRACTED_IN_MEMORY'; end
R.boldFile=boldFile;
R.qualityMaskFile=qualityFile;
R.roiFiles=roiFiles;
R.roiTypes=roiTypes;
R.TR=TR;
R.N=N;
R.timeSeconds=(0:N-1)'*TR;
R.timeseries=series;
R.coverage=coverage;
R.qualitySelection=qualitySelection;
R.outputDir='';

if cfg.writeOutputs
    out=local_path(cfg.outputDir);
    assert(~isempty(out),'outputDir required when writeOutputs=true');
    assert(~isfolder(out) && ~isfile(out), ...
        'Refusing existing output path: %s',out);
    parent=fileparts(out);
    assert(isfolder(parent),'Output parent must already exist: %s',parent);
    outCanonical=char(java.io.File(out).getCanonicalPath());
    inputFiles=[{boldFile},roiFiles];
    if ~isempty(qualityFile), inputFiles{end+1}=qualityFile; end
    for ki=1:numel(inputFiles)
        inCanonical=char(java.io.File(inputFiles{ki}).getCanonicalPath());
        assert(~strcmp(outCanonical,inCanonical) && ...
            ~startsWith(inCanonical,[outCanonical filesep]), ...
            'Output directory may not contain an input');
    end
    assert(exist('cicada_task_glm_sha256','file')==2, ...
        'SHA256 helper not on MATLAB path');
    hashes=cellfun(@cicada_task_glm_sha256,inputFiles, ...
        'UniformOutput',false);
    stage=tempname(parent);
    assert(~isfolder(stage) && ~isfile(stage),'Staging collision');
    [ok,msg]=mkdir(stage);
    assert(ok,'Cannot create staging directory: %s',msg);
    cleanStage=onCleanup(@()local_remove_stage(stage)); %#ok<NASGU>

    T=array2table([R.timeSeconds series], ...
        'VariableNames',[{'time_seconds'};labels]);
    writetable(T,fullfile(stage,'roi_mean_timeseries.csv'));
    writetable(coverage,fullfile(stage,'roi_coverage.csv'));

    manifest=struct();
    manifest.schema='cicada_roi_timeseries_v02';
    manifest.input_mode='user_supplied_rois';
    manifest.bold_file=boldFile;
    manifest.bold_sha256=hashes{1};
    manifest.quality_mask_file=qualityFile;
    manifest.quality_selection=qualitySelection;
    manifest.roi_files=roiFiles;
    manifest.roi_sha256=hashes(2:1+numel(roiFiles));
    if ~isempty(qualityFile)
        manifest.quality_mask_sha256=hashes{end};
    else
        manifest.quality_mask_sha256='';
    end
    manifest.TR_seconds=TR;
    manifest.n_volumes=N;
    manifest.n_rois=K;
    manifest.min_coverage_fraction=cfg.minCoverageFraction;
    manifest.min_valid_voxels=cfg.minValidVoxels;
    manifest.valid_voxel_rule='all_finite_all_nonzero_temporally_nonconstant';
    manifest.roi_column_names=labels;
    manifest.roi_label_ids=labelIDs;
    manifest.roi_types=roiTypes;
    manifest.roi_type_per_column=roiKinds;
    manifest.roi_volume_index=volumeIndices;
    manifest.roi_weight_definition='finite_nonnegative_renormalized_after_quality_signal';
    manifest.coverage_columns='voxels_and_weight_mass';
    manifest.weighted_min_coverage_metric='retained_weight_mass_fraction';
    manifest.empty_roi_policy='NaN_time_series_status_empty_roi';
    manifest.output_status='descriptive_means_no_inference';
    fp=fopen(fullfile(stage,'roi_provenance.json'),'w');
    assert(fp>=0,'Could not create staged provenance JSON');
    bytes=fwrite(fp,jsonencode(manifest),'char');
    closeStatus=fclose(fp);
    assert(bytes>0 && closeStatus==0,'Could not write provenance JSON');
    assert(~isfolder(out) && ~isfile(out),'Output appeared during staging');
    [ok,msg]=movefile(stage,out);
    assert(ok,'Cannot promote ROI output: %s',msg);
    R.status='EXPORTED_NO_OVERWRITE';
    R.outputDir=out;
end
fprintf('ROI %s | %d parcels, %d volumes, %d ready time series | mask=%s\n', ...
    R.status,K,N,nnz(status=="ok"),qualitySelection);
end

function [labels,files,labelIDs,volumes,kinds,indices,weights] = ...
    local_append_roi(labels,files,labelIDs,volumes,kinds,indices,weights, ...
    file,ids,w,code,volume,kind)
labels{end+1,1}=sprintf('roi_%03d',numel(labels)+1);
files{end+1,1}=file;
labelIDs(end+1,1)=code;
volumes(end+1,1)=volume;
kinds{end+1,1}=kind;
indices{end+1,1}=ids;
weights{end+1,1}=w;
end
function fraction=local_safe_fraction(numerator,denominator)
fraction=NaN(size(numerator));
positive=denominator>0;
fraction(positive)=numerator(positive)./denominator(positive);
end

function cfg=local_def(cfg,name,value)
if ~isfield(cfg,name), cfg.(name)=value; end
end
function tf=local_flag(v)
tf=isscalar(v) && (islogical(v) || ...
    (isnumeric(v) && isfinite(v) && ismember(v,[0 1])));
end
function path=local_path(v)
assert((ischar(v) && isrow(v)) || ...
    (isstring(v) && isscalar(v)), 'Expected scalar path string');
path=char(v);
assert(~isempty(strtrim(path)),'Path cannot be empty');
end
function local_same_grid(ref,subject,file,allow4D)
a=double(ref.ImageSize); b=double(subject.ImageSize);
assert((numel(b)==3 || (allow4D && numel(b)==4 && b(4)>=1)) ...
    && isequal(a(1:3),b(1:3)), ...
    'ROI/quality image dimensions differ from BOLD: %s',file);
assert(all(abs(double(ref.PixelDimensions(1:3))- ...
    double(subject.PixelDimensions(1:3)))<1e-5), ...
    'ROI/quality voxel sizes differ from BOLD: %s',file);
A=local_affine(ref); B=local_affine(subject);
assert(max(abs(A(:)-B(:)))<1e-4, ...
    'ROI/quality affine differs from BOLD: %s',file);
if isfield(ref,'SpaceUnits') && isfield(subject,'SpaceUnits')
    assert(strcmpi(char(ref.SpaceUnits),char(subject.SpaceUnits)), ...
        'Spatial units differ: %s',file);
end
end
function T=local_affine(info)
assert(isfield(info,'Transform') && ~isempty(info.Transform), ...
    'NIfTI image missing spatial transform');
t=info.Transform;
if isobject(t) && isprop(t,'T')
    T=double(t.T);
elseif isstruct(t) && isfield(t,'T')
    T=double(t.T);
else
    error('Unsupported NIfTI spatial transform');
end
assert(isequal(size(T),[4 4]) && all(isfinite(T(:))) && ...
    rcond(T)>1e-12,'Invalid NIfTI spatial transform');
end
function local_remove_stage(path)
if isfolder(path), rmdir(path,'s'); end
end
