function manifest = cicada_group_merge_support_products(product_structs, image_keep, output_dir)
% CICADA_GROUP_MERGE_SUPPORT_PRODUCTS Merge retained support products.
if isempty(product_structs)
    manifest = table();
    return
end
if ~iscell(product_structs)
    error('CICADA:SupportProductStructs', 'product_structs must be a cell array.');
end
image_keep = logical(image_keep(:));
if numel(image_keep) ~= numel(product_structs)
    error('CICADA:SupportKeepLength', 'image_keep length does not match support products.');
end
retained = product_structs(image_keep);
if isempty(retained)
    manifest = table();
    return
end

fields = {'DirectReliabilityMask','RobustTSNR','RobustTSNRNormalized','MeanNormalized', ...
    'SignalICOverlapRaw','NoiseICOverlapRaw','SignalCoreP50','SignalCoreP95', ...
    'NoiseCoreP50','NoiseCoreP95','SignalNoiseCategoryP50', ...
    'SignalNoiseBalance','AnyICEvidence','SignalICOverlapSmoothed','NoiseICOverlapSmoothed', ...
    'SignalExtendedAdaptive','SignalExtendedFixed'};
output_names = {'direct_reliability_masks','robust_tsnr_maps','robust_tsnr_normalized_maps', ...
    'mean_normalized_maps','signal_ic_overlap_raw_maps','noise_ic_overlap_raw_maps', ...
    'signal_core_p50_masks','signal_core_p95_masks','noise_core_p50_masks', ...
    'noise_core_p95_masks','signal_noise_category_p50_maps', ...
    'signal_minus_noise_evidence_maps','maximum_any_ic_evidence_maps', ...
    'signal_ic_overlap_smoothed_maps','noise_ic_overlap_smoothed_maps', ...
    'signal_extended_adaptive_masks','signal_extended_fixed_masks'};
aggregate_modes = {'both','median','median','median','median','median','fraction','fraction', ...
    'fraction','fraction','none','median','median','median','median','fraction','fraction'};

release_info = cicada_support_release_info();
manifest = table((1:numel(retained))', 'VariableNames', {'RetainedIndex'});
manifest.SupportRelease = repmat(string(release_info.ReleaseTag), numel(retained), 1);
manifest.BaselineCommit = repmat(string(release_info.BaselineCommit), numel(retained), 1);
manifest.TaskDirectory = strings(numel(retained),1);
manifest.CICADAType = strings(numel(retained),1);
for i = 1:numel(retained)
    if isstruct(retained{i})
        if isfield(retained{i}, 'TaskDirectory')
            manifest.TaskDirectory(i) = string(retained{i}.TaskDirectory);
        end
        if isfield(retained{i}, 'CICADAType')
            manifest.CICADAType(i) = string(retained{i}.CICADAType);
        end
    end
end
for f = 1:numel(fields)
    paths = strings(numel(retained),1);
    for i = 1:numel(retained)
        if isstruct(retained{i}) && isfield(retained{i}, fields{f})
            paths(i) = string(retained{i}.(fields{f}));
        end
    end
    manifest.(fields{f}) = paths;
    valid = strlength(paths) > 0 & cellfun(@isfile, cellstr(paths));
    if ~all(valid)
        fprintf('Skipping group support product %s: available for %d/%d retained images.\n', ...
            fields{f}, nnz(valid), numel(valid));
        continue
    end
    stack_path = fullfile(output_dir, [output_names{f}, '.nii.gz']);
    quoted = cellfun(@(x) ['"', x, '"'], cellstr(paths), 'UniformOutput', false);
    [status, msg] = call_fsl(sprintf('fslmerge -a "%s" %s', stack_path, strjoin(quoted, ' ')));
    if status ~= 0
        error('CICADA:FSLFailure', 'Failed to merge %s: %s', fields{f}, string(msg));
    end
    switch aggregate_modes{f}
        case 'fraction'
            out = fullfile(output_dir, ['group_', output_names{f}, '_fraction.nii.gz']);
            local_fslmaths(stack_path, '-Tmean', out);
        case 'median'
            out = fullfile(output_dir, ['group_', output_names{f}, '_median.nii.gz']);
            local_fslmaths(stack_path, '-Tmedian', out);
        case 'both'
            out1 = fullfile(output_dir, ['group_', output_names{f}, '_intersection.nii.gz']);
            out2 = fullfile(output_dir, ['group_', output_names{f}, '_fraction.nii.gz']);
            local_fslmaths(stack_path, '-Tmin', out1);
            local_fslmaths(stack_path, '-Tmean', out2);
    end
end
writetable(manifest, fullfile(output_dir, 'support_products_manifest.tsv'), ...
    'FileType', 'text', 'Delimiter', '\t');
end

function local_fslmaths(input_path, operation, output_path)
[status, msg] = call_fsl(sprintf('fslmaths "%s" %s "%s"', input_path, operation, output_path));
if status ~= 0
    error('CICADA:FSLFailure', 'Failed group aggregation: %s', string(msg));
end
end
