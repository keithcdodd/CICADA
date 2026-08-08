function manifest = cicada_group_merge_support_products(product_structs, image_keep, output_dir)
% CICADA_GROUP_MERGE_SUPPORT_PRODUCTS Merge retained support products.
%
% New Group CICADA support outputs are organized as:
%
%   support_products/
%       core/         Primary downstream support products
%       qc/           QC and sensitivity products
%       provenance/   Manifest/provenance
%
% This helper does not modify legacy Group CICADA outputs such as
% group_funcmask.nii.gz or group_signal_funcmask.nii.gz.

if isempty(product_structs)
    manifest = table();
    return
end

if ~iscell(product_structs)
    error('CICADA:SupportProductStructs', ...
        'product_structs must be a cell array.');
end

image_keep = logical(image_keep(:));

if numel(image_keep) ~= numel(product_structs)
    error('CICADA:SupportKeepLength', ...
        'image_keep length does not match support products.');
end

retained = product_structs(image_keep);

if isempty(retained)
    manifest = table();
    return
end


%% Set up organized output folders

support_root = fullfile(output_dir, 'support_products');
stage_root = fullfile(output_dir, 'support_products_staging');

% Build the complete new support-products tree in staging first.
% This avoids deleting a previously successful support_products directory
% if generation fails partway through.
if isfolder(stage_root)
    rmdir(stage_root, 's');
end

core_dir = fullfile(stage_root, 'core');
qc_dir = fullfile(stage_root, 'qc');
provenance_dir = fullfile(stage_root, 'provenance');

mkdir(core_dir);
mkdir(qc_dir);
mkdir(provenance_dir);


%% Product definitions

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

% These are the primary downstream support products.
% Everything else remains available, but is classified as QC/sensitivity.
core_fields = {'DirectReliabilityMask','SignalICOverlapRaw'};


%% Manifest

release_info = cicada_support_release_info();

manifest = table((1:numel(retained))', ...
    'VariableNames', {'RetainedIndex'});

manifest.SupportRelease = repmat(string(release_info.ReleaseTag), ...
    numel(retained), 1);

manifest.BaselineCommit = repmat(string(release_info.BaselineCommit), ...
    numel(retained), 1);

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


%% Merge support products

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

        fprintf(['Skipping group support product %s: available for ' ...
            '%d/%d retained images.\n'], ...
            fields{f}, nnz(valid), numel(valid));

        continue

    end


    % Primary products go in core/. Everything else is QC/sensitivity.

    if ismember(fields{f}, core_fields)
        target_dir = core_dir;
    else
        target_dir = qc_dir;
    end


    stack_path = fullfile(target_dir, ...
        [output_names{f}, '.nii.gz']);

    quoted = cellfun(@(x) ['"', x, '"'], ...
        cellstr(paths), 'UniformOutput', false);

    [status, msg] = call_fsl(sprintf( ...
        'fslmerge -a "%s" %s', ...
        stack_path, strjoin(quoted, ' ')));

    if status ~= 0
        error('CICADA:FSLFailure', ...
            'Failed to merge %s: %s', ...
            fields{f}, string(msg));
    end


    %% Group aggregate

    switch aggregate_modes{f}

        case 'fraction'

            out = fullfile(target_dir, ...
                ['group_', output_names{f}, '_fraction.nii.gz']);

            local_fslmaths(stack_path, '-Tmean', out);


        case 'median'

            out = fullfile(target_dir, ...
                ['group_', output_names{f}, '_median.nii.gz']);

            local_fslmaths(stack_path, '-Tmedian', out);


        case 'both'

            out1 = fullfile(target_dir, ...
                ['group_', output_names{f}, '_intersection.nii.gz']);

            out2 = fullfile(target_dir, ...
                ['group_', output_names{f}, '_fraction.nii.gz']);

            local_fslmaths(stack_path, '-Tmin', out1);
            local_fslmaths(stack_path, '-Tmean', out2);

    end

end


%% Write provenance manifest

writetable(manifest, ...
    fullfile(provenance_dir, 'support_products_manifest.tsv'), ...
    'FileType', 'text', ...
    'Delimiter', '\t');


%% Write a simple human-readable guide

local_write_output_guide(stage_root);


%% Finalize support-products directory

% Only replace the canonical directory after all products, aggregates,
% provenance, and documentation have been generated successfully.
if isfolder(support_root)
    rmdir(support_root, 's');
end

[status, msg] = movefile(stage_root, support_root);

if ~status
    error('CICADA:SupportOutputMoveFailure', ...
        'Could not finalize group support-products directory: %s', msg);
end


%% Console guidance

fprintf('\n');
fprintf('Group CICADA support products:\n');
fprintf('  Primary products:  support_products/core/\n');
fprintf('  QC/sensitivity:     support_products/qc/\n');
fprintf('  Provenance:         support_products/provenance/\n');
fprintf('  Output guide:       support_products/README_OUTPUTS.txt\n');
fprintf('\n');

end


function local_fslmaths(input_path, operation, output_path)

[status, msg] = call_fsl(sprintf( ...
    'fslmaths "%s" %s "%s"', ...
    input_path, operation, output_path));

if status ~= 0
    error('CICADA:FSLFailure', ...
        'Failed group aggregation: %s', string(msg));
end

end


function local_write_output_guide(support_root)

readme_path = fullfile(support_root, 'README_OUTPUTS.txt');

fid = fopen(readme_path, 'w');

if fid == -1
    warning('CICADA:SupportReadmeWriteFailure', ...
        'Could not write support-products output guide.');
    return
end

cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>


fprintf(fid, 'GROUP CICADA SUPPORT PRODUCTS\n');
fprintf(fid, '=============================\n\n');


fprintf(fid, 'START HERE\n');
fprintf(fid, '----------\n\n');

fprintf(fid, 'core/direct_reliability_masks.nii.gz\n');
fprintf(fid, ['  4-D binary stack, one volume per retained image. Indicates where\n' ...
    '  direct functional data passed CICADA reliability criteria.\n\n']);

fprintf(fid, 'core/group_direct_reliability_masks_fraction.nii.gz\n');
fprintf(fid, ['  3-D group map. Each voxel gives the fraction of retained images\n' ...
    '  with directly reliable functional data at that voxel.\n\n']);

fprintf(fid, 'core/group_direct_reliability_masks_intersection.nii.gz\n');
fprintf(fid, ['  3-D strict intersection. A voxel survives only if direct functional\n' ...
    '  data were reliable in every retained image.\n\n']);

fprintf(fid, 'core/signal_ic_overlap_raw_maps.nii.gz\n');
fprintf(fid, ['  4-D continuous retained-signal IC spatial-evidence stack.\n' ...
    '  This is evidence, not a calibrated probability of neuronal signal.\n\n']);

fprintf(fid, 'core/group_signal_ic_overlap_raw_maps_median.nii.gz\n');
fprintf(fid, ['  3-D median retained-signal IC evidence across retained images.\n\n']);


fprintf(fid, 'HOW TO THINK ABOUT THESE PRODUCTS\n');
fprintf(fid, '---------------------------------\n\n');

fprintf(fid, ['Direct reliability asks whether sufficiently usable functional\n' ...
    'measurements are available at a voxel.\n\n']);

fprintf(fid, ['Signal/noise IC evidence asks how strongly components classified as\n' ...
    'signal or noise spatially support that voxel.\n\n']);

fprintf(fid, ['These are complementary concepts. Signal and noise IC evidence may\n' ...
    'coexist at the same voxel. Noise evidence is not automatically an\n' ...
    'exclusion mask.\n\n']);


fprintf(fid, 'QC AND SENSITIVITY PRODUCTS\n');
fprintf(fid, '---------------------------\n\n');

fprintf(fid, ['Files in qc/ include tSNR-derived maps, noise IC evidence,\n' ...
    'thresholded signal/noise cores, smoothed evidence, extended masks,\n' ...
    'and signal-minus-noise summaries.\n\n']);

fprintf(fid, ['They are retained for QC, interpretation, and sensitivity analyses.\n' ...
    'They should not be viewed as competing default analysis masks.\n\n']);


fprintf(fid, 'PROVENANCE\n');
fprintf(fid, '----------\n\n');

fprintf(fid, ['provenance/support_products_manifest.tsv records source products\n' ...
    'and retained-image provenance for the group products.\n\n']);


fprintf(fid, 'LEGACY GROUP CICADA MASKS\n');
fprintf(fid, '-------------------------\n\n');

fprintf(fid, ['Existing Group CICADA products such as group_funcmask.nii.gz and\n' ...
    'group_signal_funcmask.nii.gz remain outside this folder for backward\n' ...
    'compatibility and are not replaced by these support products.\n']);

end