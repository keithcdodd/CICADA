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

required = {products.TemporalMean, products.TemporalStd, products.TemporalMedian, ...
    products.TemporalMAD, products.RobustTSNR, products.MeanNormalized, ...
    products.RobustTSNRNormalized, products.DirectReliabilityMask, products.Provenance};
if options.WriteStandardTSNR
    required{end+1} = products.StandardTSNR; %#ok<AGROW>
end
if ~options.Overwrite && all(cellfun(@isfile, required))
    return
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

robust_sigma = options.RobustSigmaScale .* tmad;
robust_tsnr = zeros(size(tmedian), 'single');
valid_robust = mask & isfinite(tmedian) & isfinite(robust_sigma) & robust_sigma > 0;
robust_tsnr(valid_robust) = tmedian(valid_robust) ./ robust_sigma(valid_robust);
robust_tsnr(~isfinite(robust_tsnr) | robust_tsnr < 0) = 0;
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

outlier_cutoff = 0;
if options.UseLowRobustTSNROutlierExclusion
    log_vals = log1p(valid_rtsnr_vals);
    low_flags = isoutlier(log_vals, 'median') & log_vals < median(log_vals, 'omitnan');
    if any(low_flags)
        outlier_cutoff = max(valid_rtsnr_vals(low_flags));
    end
end
fixed_cutoff = 0;
if ~isempty(options.MinimumRobustTSNR)
    fixed_cutoff = options.MinimumRobustTSNR;
end
final_cutoff = max(outlier_cutoff, fixed_cutoff);
reliable_mask = valid_robust & robust_tsnr > final_cutoff;
cicada_write_nifti_like(uint8(reliable_mask), funcmask_constrained, products.DirectReliabilityMask, 'uint8');

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
    "constrained_mask_plus_finite_nonconstant_timeseries_plus_optional_low_robust_tSNR_exclusion", ...
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

