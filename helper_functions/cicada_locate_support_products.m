function products = cicada_locate_support_products(task_dir, cicada_type)
% CICADA_LOCATE_SUPPORT_PRODUCTS Resolve current and legacy support products.
products = struct();
products.TaskDirectory = task_dir;
products.CICADAType = char(lower(string(cicada_type)));

products.FuncMaskConstrained = local_existing(fullfile(task_dir, 'funcmask_constrained.nii.gz'));
direct_dir = fullfile(task_dir, 'support_products', 'direct_data');
products.DirectReliabilityMask = local_existing(fullfile(direct_dir, 'funcmask_direct_reliability.nii.gz'));
products.RobustTSNR = local_existing(fullfile(direct_dir, 'funcdata_robust_tsnr.nii.gz'));
products.RobustTSNRNormalized = local_existing(fullfile(direct_dir, 'funcdata_robust_tsnr_normalized.nii.gz'));
products.MeanNormalized = local_existing(fullfile(direct_dir, 'funcdata_mean_normalized.nii.gz'));

ic_dir = fullfile(task_dir, 'support_products', ['ic_', products.CICADAType]);
select_dir = fullfile(task_dir, ['ic_', products.CICADAType, '_selection']);
products.SignalICOverlapRaw = local_first_existing({ ...
    fullfile(select_dir, 'SignalICOverlap.nii.gz'), fullfile(task_dir, 'SignalICOverlap.nii.gz')});
products.NoiseICOverlapRaw = local_first_existing({ ...
    fullfile(select_dir, 'NoiseICOverlap.nii.gz'), fullfile(task_dir, 'NoiseICOverlap.nii.gz')});
products.SignalCoreP50 = local_existing(fullfile(ic_dir, 'signal_ic_evidence_core_p50.nii.gz'));
products.SignalCoreP95 = local_existing(fullfile(ic_dir, 'signal_ic_evidence_core_p95.nii.gz'));
products.NoiseCoreP50 = local_existing(fullfile(ic_dir, 'noise_ic_evidence_core_p50.nii.gz'));
products.NoiseCoreP95 = local_existing(fullfile(ic_dir, 'noise_ic_evidence_core_p95.nii.gz'));
products.SignalNoiseCategoryP50 = local_existing(fullfile(ic_dir, 'signal_noise_evidence_category_p50.nii.gz'));
products.SignalNoiseBalance = local_existing(fullfile(ic_dir, 'signal_minus_noise_ic_evidence.nii.gz'));
products.AnyICEvidence = local_existing(fullfile(ic_dir, 'maximum_signal_or_noise_ic_evidence.nii.gz'));
products.SignalICOverlapSmoothed = local_pattern_existing(ic_dir, 'signal_ic_overlap_fwhm*.nii.gz');
products.NoiseICOverlapSmoothed = local_pattern_existing(ic_dir, 'noise_ic_overlap_fwhm*.nii.gz');
products.SignalExtendedAdaptive = local_pattern_existing(ic_dir, 'signal_ic_evidence_extended_fwhm*_thrP*.nii.gz');
products.SignalExtendedFixed = local_pattern_existing(ic_dir, 'signal_ic_evidence_extended_fwhm*_p*.nii.gz');
products.LegacySignalMask = local_first_existing({ ...
    fullfile(task_dir, sprintf('funcmask_CICADA_%s_signal_constrained.nii.gz', products.CICADAType)), ...
    fullfile(select_dir, sprintf('funcmask_CICADA_%s_signal_constrained.nii.gz', products.CICADAType))});
end

function path = local_existing(candidate)
if isfile(candidate)
    path = candidate;
else
    path = '';
end
end

function path = local_first_existing(candidates)
path = '';
for i = 1:numel(candidates)
    if isfile(candidates{i})
        path = candidates{i};
        return
    end
end
end

function path = local_pattern_existing(folder, pattern)
path = '';

if ~isfolder(folder)
    return
end

hits = dir(fullfile(folder, pattern));

if isscalar(hits)

    path = fullfile(hits(1).folder, hits(1).name);

elseif numel(hits) > 1

    % Multiple sensitivity/QC derivatives may exist if support-product
    % settings have changed across runs. Do not silently choose one,
    % because that could mix different parameterizations across subjects.
    hit_names = sort({hits.name});
    match_text = strjoin(hit_names, '\n  ');

    warning('CICADA:AmbiguousSupportProduct', ...
        ['Multiple support products match pattern "%s" in:\n' ...
         '  %s\n' ...
         'Skipping this optional support product rather than choosing ' ...
         'one silently.\n' ...
         'Matches:\n' ...
         '  %s'], ...
        pattern, folder, match_text);

end

end
