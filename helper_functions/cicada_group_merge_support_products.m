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

% Step 9: suppress a Group direct-mask duplicate only after verifying the
% retained per-scan masks AND the ordered canonical Group 4-D stack.
% Failure to prove equality takes the legacy path (write direct products).
omit_duplicate_direct = local_group_direct_masks_duplicate(retained, output_dir);
direct_group_status = "not_written_missing_sources";


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

    if strcmp(fields{f}, 'DirectReliabilityMask') && omit_duplicate_direct
        direct_group_status = "omitted_verified_duplicate";
        fprintf(['Skipping the three redundant Group direct-reliability ', ...
            'outputs: canonical Group masks verified identical.\n']);
        continue
    end

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

    if strcmp(fields{f}, 'DirectReliabilityMask')
        direct_group_status = "written";
    end

end

% This column describes group outputs; the DirectReliabilityMask column
% above always records retained per-scan source paths, even if omitted.
manifest.GroupDirectReliabilityOutputs = repmat( ...
    direct_group_status, numel(retained), 1);

%% Write provenance manifest

writetable(manifest, ...
    fullfile(provenance_dir, 'support_products_manifest.tsv'), ...
    'FileType', 'text', ...
    'Delimiter', '\t');


%% Write a simple human-readable guide

local_write_output_guide(stage_root, direct_group_status);


%% Finalize support-products directory

% Regenerate by replacement, with restoration of the previous successful
% tree on an ordinary promotion failure. This does not alter output contents.
cicada_group_finalize_support_products(stage_root, support_root);


%% Console guidance

fprintf('\n');
fprintf('Group CICADA support products:\n');
fprintf('  Primary products:  support_products/core/\n');
fprintf('  QC/sensitivity:     support_products/qc/\n');
fprintf('  Provenance:         support_products/provenance/\n');
fprintf('  Output guide:       support_products/README_OUTPUTS.txt\n');
fprintf('\n');

end


function yes = local_group_direct_masks_duplicate(retained, output_dir)
% Confirm exact binary equivalence in manifest/Group 4-D order. Use fslroi
% to read one Group volume at a time; never load the entire 4-D stack.
yes = false;
stack_path = fullfile(output_dir, 'funcmasks.nii.gz');
if ~isfile(stack_path)
    return
end

try
    info = niftiinfo(stack_path);
    n_volumes = 1;
    if numel(info.ImageSize) >= 4
        n_volumes = double(info.ImageSize(4));
    end
    strict_path = fullfile(output_dir, 'group_funcmask.nii.gz');
    fraction_path = fullfile(output_dir, 'group_funcmask_fraction.nii.gz');
    if n_volumes ~= numel(retained) || ...
            ~isfile(strict_path) || ~isfile(fraction_path)
        return
    end
    expected_strict = [];
    expected_fraction_sum = [];

    % Keep temporary extraction outside all CICADA output trees.
    temporary_dir = tempname;
    mkdir(temporary_dir);
    temporary_cleanup = onCleanup(@() local_remove_temp_dir(temporary_dir)); %#ok<NASGU>

    for i = 1:numel(retained)
        p = retained{i};
        if ~isstruct(p) || ~isfield(p, 'TaskDirectory') || ...
                ~isfield(p, 'FuncMaskConstrained') || ...
                ~isfield(p, 'DirectReliabilityMask')
            return
        end
        task_dir = char(string(p.TaskDirectory));
        constrained_path = char(string(p.FuncMaskConstrained));
        direct_path = char(string(p.DirectReliabilityMask));
        if isempty(task_dir) || isempty(constrained_path) || ...
                isempty(direct_path) || ~isfile(constrained_path) || ...
                ~isfile(direct_path)
            return
        end
        % Unified provenance is required; its voxel count must agree.
        provenance_path = fullfile(task_dir, ...
            'funcmask_constrained_provenance.tsv');
        if ~isfile(provenance_path)
            return
        end
        P = readtable(provenance_path, 'FileType', 'text', ...
            'Delimiter', '\t');
        if height(P) ~= 1 || ...
                ~ismember('FinalVoxels', P.Properties.VariableNames)
            return
        end

        % Extract the actual canonical volume. Its index must correspond to
        % RetainedIndex=i, not the original pre-exclusion scan number.
        extracted_path = fullfile(temporary_dir, ...
            sprintf('canonical_%04d.nii.gz', i));
        [status, ~] = call_fsl(sprintf('fslroi "%s" "%s" %d 1', ...
            stack_path, extracted_path, i-1));
        if status ~= 0 || ~isfile(extracted_path)
            return
        end
        cicada_assert_nifti_geometry(constrained_path, ...
            direct_path, stack_path, extracted_path, ...
            strict_path, fraction_path);
        constrained = double(niftiread(constrained_path));
        direct = double(niftiread(direct_path));
        group_volume = double(niftiread(extracted_path));
        if ~isequal(size(constrained), size(direct), size(group_volume)) || ...
                ~all(isfinite(constrained(:))) || ...
                ~all(isfinite(direct(:))) || ...
                ~all(isfinite(group_volume(:))) || ...
                ~all(constrained(:)==0 | constrained(:)==1) || ...
                ~all(direct(:)==0 | direct(:)==1) || ...
                ~all(group_volume(:)==0 | group_volume(:)==1) || ...
                ~isequal(constrained, direct, group_volume) || ...
                ~isnumeric(P.FinalVoxels) || ...
                P.FinalVoxels(1) ~= nnz(constrained)
            return
        end
        if i == 1
            expected_strict = logical(constrained);
            expected_fraction_sum = constrained;
        else
            expected_strict = expected_strict & logical(constrained);
            expected_fraction_sum = expected_fraction_sum + constrained;
        end
    end

    % Check that the retained 3-D canonical summaries actually match the
    % same ordered stack, not just that files with those names exist.
    actual_strict = double(niftiread(strict_path));
    actual_fraction = double(niftiread(fraction_path));
    if ~isequal(size(actual_strict), size(expected_strict)) || ...
            ~isequal(size(actual_fraction), size(expected_strict)) || ...
            ~all(isfinite(actual_strict(:))) || ...
            ~all(isfinite(actual_fraction(:))) || ...
            ~isequal(actual_strict, double(expected_strict)) || ...
            any(abs(actual_fraction(:) - ...
                expected_fraction_sum(:) / numel(retained)) > 1e-6)
        return
    end
    yes = true;
catch ME
    % A verification error must not silently drop a potentially distinct
    % legacy product. The old Group-direct output remains the default.
    warning('CICADA:GroupDirectMaskVerification', ...
        'Cannot prove duplicate Group masks; retaining legacy outputs: %s', ...
        ME.message);
end
end

function local_remove_temp_dir(path)
if isfolder(path)
    rmdir(path, 's');
end
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


function local_write_output_guide(support_root, direct_group_status)

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

if direct_group_status == "omitted_verified_duplicate"
    fprintf(fid, 'VERIFIED EQUIVALENT GROUP MASKS (no duplicate direct products)\n');
    fprintf(fid, ['  funcmasks.nii.gz (4-D retained-order stack)\n' ...
        '  group_funcmask.nii.gz (3-D strict intersection)\n' ...
        '  group_funcmask_fraction.nii.gz (3-D retained-scan fraction)\n' ...
        '  These canonical files are in the parent Group directory.\n' ...
        '  They were checked against every retained direct mask.\n\n']);
elseif direct_group_status == "written"
    fprintf(fid, 'core/direct_reliability_masks.nii.gz\n');
    fprintf(fid, ['  4-D binary stack, one volume per retained image. Indicates where\n' ...
        '  direct functional data passed CICADA reliability criteria.\n\n']);
    fprintf(fid, 'core/group_direct_reliability_masks_fraction.nii.gz\n');
    fprintf(fid, ['  3-D group map. Each voxel gives the fraction of retained images\n' ...
        '  with directly reliable functional data at that voxel.\n\n']);
    fprintf(fid, 'core/group_direct_reliability_masks_intersection.nii.gz\n');
    fprintf(fid, ['  3-D strict intersection. A voxel survives only if direct functional\n' ...
        '  data were reliable in every retained image.\n\n']);
else
    fprintf(fid, ['Group direct-reliability maps were not written because not all\n' ...
        '  retained direct-mask source files were available.\n\n']);
end

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