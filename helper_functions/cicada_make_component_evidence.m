function products = cicada_make_component_evidence(task_dir, cicada_type, signal_overlap_path, noise_overlap_path, support_options)
% CICADA_MAKE_COMPONENT_EVIDENCE Create explicit ICA-derived evidence products.
%
% Raw SignalICOverlap/NoiseICOverlap remain the maximum MELODIC probability
% across ICs labeled signal/noise. This helper adds fixed-threshold local
% cores and a modestly smoothed regional-evidence map without replacing the
% historical sigma=6 mm / -thrP 50 signal-constrained mask.

if nargin < 5
    support_options = struct();
end
options = cicada_support_options(support_options);
products = struct();
if ~options.Enabled
    return
end

cicada_type = lower(string(cicada_type));
if ~ismember(cicada_type, ["auto","manual"])
    error('CICADA:InvalidCICADAType', 'cicada_type must be auto or manual.');
end
if ~isfile(signal_overlap_path)
    error('CICADA:MissingSignalICOverlap', 'Signal IC overlap file does not exist: %s', signal_overlap_path);
end
if ~isfile(noise_overlap_path)
    error('CICADA:MissingNoiseICOverlap', 'Noise IC overlap file does not exist: %s', noise_overlap_path);
end

constrained_mask = fullfile(task_dir, 'funcmask_constrained.nii.gz');
if ~isfile(constrained_mask)
    fallback_mask = fullfile(task_dir, 'funcmask.nii.gz');
    if ~isfile(fallback_mask)
        error('CICADA:MissingFunctionalMask', 'No constrained or original functional mask found.');
    end
    constrained_mask = fallback_mask;
end

out_dir = fullfile(task_dir, 'support_products', ['ic_', char(cicada_type)]);
if ~isfolder(out_dir)
    mkdir(out_dir);
end
products.OutputDirectory = out_dir;
products.SignalICOverlapRaw = signal_overlap_path;
products.NoiseICOverlapRaw = noise_overlap_path;
products.ConstrainedMask = constrained_mask;

expected = {fullfile(out_dir, 'signal_minus_noise_ic_evidence.nii.gz'), ...
    fullfile(out_dir, 'maximum_signal_or_noise_ic_evidence.nii.gz'), ...
    fullfile(out_dir, 'signal_noise_evidence_category_p50.nii.gz'), ...
    fullfile(out_dir, sprintf('signal_ic_overlap_fwhm%g.nii.gz', options.ExtendedSmoothingFWHMmm)), ...
    fullfile(out_dir, sprintf('noise_ic_overlap_fwhm%g.nii.gz', options.ExtendedSmoothingFWHMmm)), ...
    fullfile(out_dir, sprintf('signal_ic_evidence_extended_fwhm%g_thrP%g.nii.gz', ...
        options.ExtendedSmoothingFWHMmm, options.ExtendedThresholdPercentRobustRange)), ...
    fullfile(out_dir, 'component_evidence_provenance.tsv')};
for core_probability = options.SignalCoreProbabilities
    expected{end+1} = fullfile(out_dir, ...
        ['signal_ic_evidence_core_p', local_probability_label(core_probability), '.nii.gz']); %#ok<AGROW>
end
for core_probability = options.NoiseCoreProbabilities
    expected{end+1} = fullfile(out_dir, ...
        ['noise_ic_evidence_core_p', local_probability_label(core_probability), '.nii.gz']); %#ok<AGROW>
end
if options.WriteFixedThresholdExtendedMask
    expected{end+1} = fullfile(out_dir, sprintf( ...
        'signal_ic_evidence_extended_fwhm%g_p%s.nii.gz', ...
        options.ExtendedSmoothingFWHMmm, ...
        local_probability_label(options.FixedThresholdExtendedProbability))); %#ok<AGROW>
end
if ~options.Overwrite && all(cellfun(@isfile, expected))
    products = cicada_locate_support_products(task_dir, cicada_type);
    return
end

cicada_assert_nifti_geometry(constrained_mask, signal_overlap_path, noise_overlap_path);
signal = single(niftiread(signal_overlap_path));
noise = single(niftiread(noise_overlap_path));
mask = niftiread(constrained_mask) > 0;
if ~isequal(size(signal), size(noise), size(mask))
    error('CICADA:ComponentEvidenceGeometry', 'Signal, noise, and constrained-mask arrays do not have identical dimensions.');
end
if any(~isfinite(signal(:))) || any(~isfinite(noise(:)))
    error('CICADA:NonfiniteComponentProbability', 'Signal/noise overlap maps contain nonfinite values.');
end
if any(signal(:) < -1e-6 | signal(:) > 1+1e-6) || any(noise(:) < -1e-6 | noise(:) > 1+1e-6)
    error('CICADA:ComponentProbabilityRange', 'Signal/noise overlap values must lie within [0,1].');
end
signal = min(max(signal,0),1);
noise = min(max(noise,0),1);

% Continuous descriptive maps retain coexistence rather than using noise
% evidence as an exclusion criterion.
products.SignalNoiseBalance = fullfile(out_dir, 'signal_minus_noise_ic_evidence.nii.gz');
products.AnyICEvidence = fullfile(out_dir, 'maximum_signal_or_noise_ic_evidence.nii.gz');
balance = zeros(size(signal), 'single');
any_evidence = zeros(size(signal), 'single');
balance(mask) = signal(mask) - noise(mask);
any_evidence(mask) = max(signal(mask), noise(mask));
cicada_write_nifti_like(balance, constrained_mask, products.SignalNoiseBalance, 'single');
cicada_write_nifti_like(any_evidence, constrained_mask, products.AnyICEvidence, 'single');

for p = options.SignalCoreProbabilities
    label = local_probability_label(p);
    field = matlab.lang.makeValidName(['SignalCoreP', label]);
    path = fullfile(out_dir, ['signal_ic_evidence_core_p', label, '.nii.gz']);
    cicada_write_nifti_like(uint8(mask & signal >= p), constrained_mask, path, 'uint8');
    products.(field) = path;
end
for p = options.NoiseCoreProbabilities
    label = local_probability_label(p);
    field = matlab.lang.makeValidName(['NoiseCoreP', label]);
    path = fullfile(out_dir, ['noise_ic_evidence_core_p', label, '.nii.gz']);
    cicada_write_nifti_like(uint8(mask & noise >= p), constrained_mask, path, 'uint8');
    products.(field) = path;
end

% Four-level category at p=.50 when available: 0 neither, 1 signal only,
% 2 noise only, 3 mixed signal and noise evidence.
if any(abs(options.SignalCoreProbabilities-0.5) < 1e-12) && ...
        any(abs(options.NoiseCoreProbabilities-0.5) < 1e-12)
    s50 = mask & signal >= 0.5;
    n50 = mask & noise >= 0.5;
    category = uint8(s50) + uint8(n50).*2;
    category_path = fullfile(out_dir, 'signal_noise_evidence_category_p50.nii.gz');
    cicada_write_nifti_like(category, constrained_mask, category_path, 'uint8');
    products.SignalNoiseCategoryP50 = category_path;
end

sigma_mm = options.ExtendedSmoothingFWHMmm / 2.354820045;
products.SignalICOverlapSmoothed = fullfile(out_dir, ...
    sprintf('signal_ic_overlap_fwhm%g.nii.gz', options.ExtendedSmoothingFWHMmm));
products.NoiseICOverlapSmoothed = fullfile(out_dir, ...
    sprintf('noise_ic_overlap_fwhm%g.nii.gz', options.ExtendedSmoothingFWHMmm));
q = @(p) ['"', char(p), '"'];
[status, msg] = call_fsl(sprintf('fslmaths %s -s %.8g -mas %s %s', ...
    q(signal_overlap_path), sigma_mm, q(constrained_mask), q(products.SignalICOverlapSmoothed)));
local_require_fsl(status, msg, 'smoothed signal-IC evidence');
[status, msg] = call_fsl(sprintf('fslmaths %s -s %.8g -mas %s %s', ...
    q(noise_overlap_path), sigma_mm, q(constrained_mask), q(products.NoiseICOverlapSmoothed)));
local_require_fsl(status, msg, 'smoothed noise-IC evidence');

products.SignalExtendedAdaptive = fullfile(out_dir, sprintf( ...
    'signal_ic_evidence_extended_fwhm%g_thrP%g.nii.gz', ...
    options.ExtendedSmoothingFWHMmm, options.ExtendedThresholdPercentRobustRange));
[status, msg] = call_fsl(sprintf('fslmaths %s -thrP %.8g -bin -mas %s %s', ...
    q(products.SignalICOverlapSmoothed), options.ExtendedThresholdPercentRobustRange, ...
    q(constrained_mask), q(products.SignalExtendedAdaptive)));
local_require_fsl(status, msg, 'adaptive extended signal-IC evidence mask');

if options.WriteFixedThresholdExtendedMask
    label = local_probability_label(options.FixedThresholdExtendedProbability);
    products.SignalExtendedFixed = fullfile(out_dir, sprintf( ...
        'signal_ic_evidence_extended_fwhm%g_p%s.nii.gz', ...
        options.ExtendedSmoothingFWHMmm, label));
    smoothed_signal = single(niftiread(products.SignalICOverlapSmoothed));
    fixed_mask = mask & smoothed_signal >= options.FixedThresholdExtendedProbability;
    cicada_write_nifti_like(uint8(fixed_mask), constrained_mask, products.SignalExtendedFixed, 'uint8');
end

legacy_path = fullfile(task_dir, sprintf('funcmask_CICADA_%s_signal_constrained.nii.gz', cicada_type));
products.LegacySignalMask = legacy_path;
products.Provenance = fullfile(out_dir, 'component_evidence_provenance.tsv');
release_info = cicada_support_release_info();
prov = table(string(release_info.ReleaseTag), ...
    string(release_info.BaselineRepository), string(release_info.BaselineCommit), ...
    cicada_type, string(signal_overlap_path), string(noise_overlap_path), ...
    string(constrained_mask), options.ExtendedSmoothingFWHMmm, sigma_mm, ...
    options.ExtendedThresholdPercentRobustRange, ...
    options.FixedThresholdExtendedProbability, options.LegacySmoothingSigmaMM, ...
    options.LegacyThresholdPercentRobustRange, string(legacy_path), ...
    string(strjoin(compose('%.2f', options.SignalCoreProbabilities), ',')), ...
    string(strjoin(compose('%.2f', options.NoiseCoreProbabilities), ',')), ...
    "voxelwise_maximum_across_ICs_with_the_given_label", ...
    "signal_minus_noise; descriptive_only; not_an_exclusion_rule", true, ...
    'VariableNames', {'SupportRelease','BaselineRepository','BaselineCommit', ...
    'CICADAType','SignalICOverlapSource','NoiseICOverlapSource', ...
    'ConstrainedMask','ExtendedSmoothingFWHMmm','ExtendedSmoothingSigmaMM', ...
    'ExtendedAdaptiveThresholdPercentRobustRange','ExtendedFixedThreshold', ...
    'LegacySmoothingSigmaMM','LegacyThresholdPercentRobustRange','LegacySignalMask', ...
    'SignalCoreProbabilities','NoiseCoreProbabilities','RawICAggregation', ...
    'SignalNoiseBalanceSemantics','LegacyOutputPreserved'});
writetable(prov, products.Provenance, 'FileType', 'text', 'Delimiter', '\t');
end

function label = local_probability_label(p)
label = sprintf('%02d', round(p*100));
end

function local_require_fsl(status, msg, operation)
if status ~= 0
    error('CICADA:FSLFailure', 'FSL failed during %s: %s', operation, string(msg));
end
end
