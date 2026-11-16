function products = make_funcdata_reliability(task_dir, funcfile, funcmask_constrained, support_options)
% MAKE_FUNCDATA_RELIABILITY Create direct functional-data quality products.
%
% These products are intentionally separate from ICA-derived signal evidence.
% The binary direct-reliability mask is the constrained functional mask after
% removal of invalid/constant time series and, by default, only extreme low
% robust-tSNR outliers within that subject. No universal tSNR threshold is
% imposed unless explicitly supplied in support_options.MinimumRobustTSNR.

if nargin < 4
    support_options = struct();
end
options = cicada_support_options(support_options);
products = struct();
if ~options.Enabled
    return
end

if ~isfile(funcfile)
    error('CICADA:MissingFunctionalFile', 'Functional file does not exist: %s', funcfile);
end
if ~isfile(funcmask_constrained)
    error('CICADA:MissingConstrainedMask', 'Constrained mask does not exist: %s', funcmask_constrained);
end
if ~isfolder(task_dir)
    error('CICADA:MissingTaskDirectory', 'Task directory does not exist: %s', task_dir);
end

out_dir = fullfile(task_dir, 'support_products', 'direct_data');
if ~isfolder(out_dir)
    mkdir(out_dir);
end

products.OutputDirectory = out_dir;
products.TemporalMean = fullfile(out_dir, 'funcdata_tmean.nii.gz');
products.TemporalStd = fullfile(out_dir, 'funcdata_tstd.nii.gz');
products.TemporalMedian = fullfile(out_dir, 'funcdata_tmedian.nii.gz');
products.TemporalMAD = fullfile(out_dir, 'funcdata_temporal_mad.nii.gz');
products.StandardTSNR = fullfile(out_dir, 'funcdata_tsnr.nii.gz');
products.RobustTSNR = fullfile(out_dir, 'funcdata_robust_tsnr.nii.gz');
products.MeanNormalized = fullfile(out_dir, 'funcdata_mean_normalized.nii.gz');
products.RobustTSNRNormalized = fullfile(out_dir, 'funcdata_robust_tsnr_normalized.nii.gz');
products.DirectReliabilityMask = fullfile(out_dir, 'funcmask_direct_reliability.nii.gz');
products.Provenance = fullfile(out_dir, 'funcdata_reliability_provenance.tsv');

% Unified constrained masks already include direct-data reliability screening.
% Resolve the canonical alias and check for conflicting old duplicates BEFORE
% running FSL or writing any support product. Never delete pre-existing files.
hasUnifiedMask = isfile(fullfile(task_dir, 'funcmask_constrained_provenance.tsv'));
if hasUnifiedMask
    resolved = cicada_locate_support_products(task_dir, 'auto');
    cicada_assert_nifti_geometry(funcmask_constrained, resolved.FuncMaskConstrained);
    if ~isequal(niftiread(funcmask_constrained) > 0, ...
            niftiread(resolved.FuncMaskConstrained) > 0)
        error('CICADA:ConstrainedMaskInputMismatch', ...
            'Supplied constrained mask differs from task-folder canonical mask.');
    end
    products.DirectReliabilityMask = resolved.DirectReliabilityMask;
end

% These are outputs produced by THIS function, not the canonical unified
% mask (which is an upstream input and must never be rewritten here).
required = {products.TemporalMean, products.TemporalStd, products.TemporalMedian, ...
    products.TemporalMAD, products.RobustTSNR, products.MeanNormalized, ...
    products.RobustTSNRNormalized, products.Provenance};
if ~hasUnifiedMask
    required{end+1} = products.DirectReliabilityMask; %#ok<AGROW>
end
if options.WriteStandardTSNR
    required{end+1} = products.StandardTSNR; %#ok<AGROW>
end

% Never overwrite only some existing support products when Overwrite=false.
% Reuse a complete set after inexpensive provenance/geometry checks;
% refuse a partial set rather than silently mixing different generations.
% This cannot prove source *content* identity (legacy provenance has no hashes).
if ~options.Overwrite
    existsNow = cellfun(@isfile, required);
    if any(existsNow)
        if ~all(existsNow)
            error('CICADA:PartialReliabilityProducts', ...
                ['Some direct-data support products already exist in %s, ' ...
                 'but the requested set is incomplete. No files were rewritten. ' ...
                 'Inspect the products before explicitly allowing overwrite.'], out_dir);
        end
        P = readtable(products.Provenance, 'FileType','text', 'Delimiter','\t');
        provenanceFields = {'FunctionalFile','ConstrainedMask','MADScale', ...
            'UsedLowOutlierExclusion'};
        if height(P) == 1 && all(ismember(provenanceFields, ...
                P.Properties.VariableNames))
            savedScale = str2double(string(P.MADScale(1)));
            savedOutlier = lower(string(P.UsedLowOutlierExclusion(1)));
            savedOutlierValid = ismember(savedOutlier, ...
                ["true","false","1","0"]);
            savedOutlierOn = ismember(savedOutlier,["true","1"]);
        else
            savedScale = NaN;
            savedOutlierValid = false;
            savedOutlierOn = false;
        end
        if height(P) ~= 1 || ~all(ismember(provenanceFields, ...
                P.Properties.VariableNames)) || ...
                ~strcmp(string(P.FunctionalFile(1)),string(funcfile)) || ...
                ~strcmp(string(P.ConstrainedMask(1)),string(funcmask_constrained)) || ...
                ~isequal(savedScale,double(options.RobustSigmaScale)) || ...
                ~savedOutlierValid || ...
                ~isequal(savedOutlierOn, ...
                    logical(options.UseLowRobustTSNROutlierExclusion))
            error('CICADA:ReliabilityProvenanceMismatch', ...
                ['Existing direct-data support products have incompatible ' ...
                 'source paths or reliability settings; no files were rewritten.']);
        end
        geometryProducts = required(~strcmp(required,products.Provenance));
        cicada_assert_nifti_geometry(funcmask_constrained, geometryProducts{:});
        return
    end
end

q = @(p) ['"', char(p), '"'];
[status, msg] = call_fsl(sprintf('fslmaths %s -Tmean %s', q(funcfile), q(products.TemporalMean)));
local_require_fsl(status, msg, 'temporal mean');
[status, msg] = call_fsl(sprintf('fslmaths %s -Tstd %s', q(funcfile), q(products.TemporalStd)));
local_require_fsl(status, msg, 'temporal standard deviation');
[status, msg] = call_fsl(sprintf('fslmaths %s -Tmedian %s', q(funcfile), q(products.TemporalMedian)));
local_require_fsl(status, msg, 'temporal median');
[status, msg] = call_fsl(sprintf('fslmaths %s -sub %s -abs -Tmedian %s', ...
    q(funcfile), q(products.TemporalMedian), q(products.TemporalMAD)));
local_require_fsl(status, msg, 'temporal MAD');

mask = niftiread(funcmask_constrained) > 0;
tmean = single(niftiread(products.TemporalMean));
tstd = single(niftiread(products.TemporalStd));
tmedian = single(niftiread(products.TemporalMedian));
tmad = single(niftiread(products.TemporalMAD));

cicada_assert_nifti_geometry(funcmask_constrained, products.TemporalMean, ...
    products.TemporalStd, products.TemporalMedian, products.TemporalMAD);

standard_tsnr = zeros(size(tmean), 'single');
valid_std = mask & isfinite(tmean) & isfinite(tstd) & tstd > 0;
standard_tsnr(valid_std) = tmean(valid_std) ./ tstd(valid_std);
standard_tsnr(~isfinite(standard_tsnr) | standard_tsnr < 0) = 0;
if options.WriteStandardTSNR
    cicada_write_nifti_like(standard_tsnr, funcmask_constrained, products.StandardTSNR, 'single');
end

[reliable_mask, robust_tsnr, valid_robust, ...
    outlier_cutoff, fixed_cutoff, final_cutoff] = ...
    cicada_robust_tsnr_selection( ...
        mask, tmedian, tmad, options);
cicada_write_nifti_like(robust_tsnr, funcmask_constrained, products.RobustTSNR, 'single');

valid_mean_vals = double(tmean(mask & isfinite(tmean) & tmean > 0));
if isempty(valid_mean_vals)
    error('CICADA:NoValidMeanIntensity', 'No valid positive temporal-mean voxels in constrained mask.');
end
median_mean = median(valid_mean_vals, 'omitnan');
mean_normalized = zeros(size(tmean), 'single');
mean_normalized(mask) = tmean(mask) ./ max(single(median_mean), eps('single'));
mean_normalized(~isfinite(mean_normalized) | mean_normalized < 0) = 0;
cicada_write_nifti_like(mean_normalized, funcmask_constrained, products.MeanNormalized, 'single');

valid_rtsnr_vals = double(robust_tsnr(valid_robust & robust_tsnr > 0));
if isempty(valid_rtsnr_vals)
    error('CICADA:NoValidRobustTSNR', 'No valid robust-tSNR voxels in constrained mask.');
end
median_rtsnr = median(valid_rtsnr_vals, 'omitnan');
rtsnr_normalized = zeros(size(robust_tsnr), 'single');
rtsnr_normalized(mask) = robust_tsnr(mask) ./ max(single(median_rtsnr), eps('single'));
rtsnr_normalized(~isfinite(rtsnr_normalized) | rtsnr_normalized < 0) = 0;
cicada_write_nifti_like(rtsnr_normalized, funcmask_constrained, products.RobustTSNRNormalized, 'single');

%% Unified constrained-mask compatibility

unifiedProvFile = fullfile( ...
    task_dir, 'funcmask_constrained_provenance.tsv');

if hasUnifiedMask

    % The constrained mask has already undergone
    % mean-intensity, temporal, and robust-tSNR screening.
    upstream = readtable(unifiedProvFile, ...
        'FileType', 'text', 'Delimiter', '\t');

    assert(height(upstream) == 1 && ...
        upstream.FinalVoxels(1) == nnz(mask), ...
        'Constrained-mask provenance does not match the mask.');

    % Every retained voxel must have valid robust-tSNR data.
    assert(~any(mask(:) & ~valid_robust(:)) && ...
           ~any(mask(:) & robust_tsnr(:) <= 0), ...
        'Unified constrained mask contains invalid robust-tSNR voxels.');

    % Do not apply another adaptive exclusion.
    reliable_mask = mask;

    % Report the cutoff actually applied upstream.
    outlier_cutoff = upstream.LowRobustTSNRCutoff(1);
    fixed_cutoff = upstream.FixedRobustTSNRCutoff(1);
    final_cutoff = upstream.AppliedRobustTSNRCutoff(1);

    mask_semantics = ...
        "identical_to_unified_constrained_mask_screened_upstream";

else

    % Backward compatibility with older constrained masks.
    mask_semantics = ...
        "constrained_mask_plus_finite_nonconstant_timeseries_plus_optional_low_robust_tSNR_exclusion";

end

% In unified runs, DirectReliabilityMask resolves to the canonical file.
% Do not overwrite that source NIfTI or write a redundant physical duplicate.
% Older distinct masks continue to be written as before.
if ~hasUnifiedMask
    cicada_write_nifti_like(uint8(reliable_mask), funcmask_constrained, ...
        products.DirectReliabilityMask, 'uint8');
end

n_constrained = nnz(mask);
n_valid = nnz(valid_robust);
n_reliable = nnz(reliable_mask);
release_info = cicada_support_release_info();
provenance = table(string(release_info.ReleaseTag), ...
    string(release_info.BaselineRepository), string(release_info.BaselineCommit), ...
    string(funcfile), string(funcmask_constrained), ...
    n_constrained, n_valid, n_reliable, n_reliable/max(n_constrained,1), ...
    median_mean, median_rtsnr, outlier_cutoff, fixed_cutoff, final_cutoff, ...
    options.RobustSigmaScale, options.UseLowRobustTSNROutlierExclusion, ...
    "temporal_median/(MADScale*temporal_MAD)", ...
    mask_semantics, ...
    'VariableNames', {'SupportRelease','BaselineRepository','BaselineCommit', ...
    'FunctionalFile','ConstrainedMask','ConstrainedVoxels', ...
    'ValidRobustTSNRVoxels','ReliableVoxels','ReliableFractionOfConstrained', ...
    'MedianTemporalMean','MedianRobustTSNR','LowOutlierCutoff', ...
    'FixedMinimumRobustTSNR','AppliedRobustTSNRCutoff','MADScale', ...
    'UsedLowOutlierExclusion','RobustTSNRDefinition','DirectReliabilityMaskSemantics'});
writetable(provenance, products.Provenance, 'FileType', 'text', 'Delimiter', '\t');
end

function local_require_fsl(status, msg, operation)
if status ~= 0
    error('CICADA:FSLFailure', 'FSL failed during %s: %s', operation, string(msg));
end
end

