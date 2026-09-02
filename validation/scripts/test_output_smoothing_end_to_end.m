function test_output_smoothing_end_to_end
%TEST_OUTPUT_SMOOTHING_END_TO_END Exercise the real NIfTI smoothing path.
%   Uses only synthetic impulse data. No participant data are read or
%   written. Requires MATLAB's Image Processing Toolbox.

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(fullfile(repo_root, 'helper_functions'));

test_root = tempname;
mkdir(test_root);
cleanup_test_root = onCleanup(@() rmdir(test_root, 's')); %#ok<NASGU>

local_run_case(test_root, 'isotropic', [2 2 2], 6, 's6', true);
local_run_case(test_root, 'anisotropic', [2 2 3], 6, 's6', false);
local_run_case(test_root, 'noninteger', [2.4 2.4 2.7], 5.5, 's5p5', false);
local_run_case(test_root, 'no_smoothing', [2 2 3], 0, '', false);

fprintf('PASS: end-to-end synthetic NIfTI output-smoothing tests.\n');
end

function local_run_case(test_root, case_name, voxel_sizes_mm, ...
    fwhm_mm, expected_tag, compare_legacy)

case_dir = fullfile(test_root, case_name);
output_dir = fullfile(case_dir, 'output');
mkdir(case_dir);
mkdir(output_dir);

volume_size = [31 31 31];
n_timepoints = 4;
impulse = zeros(volume_size, 'single');
center = ceil((volume_size + 1) / 2);
impulse(center(1), center(2), center(3)) = 1;
input_data = repmat(impulse, [1 1 1 n_timepoints]);

input_base = fullfile(case_dir, 'synthetic_input.nii');
niftiwrite(input_data, input_base, 'Compressed', true);
input_file = [input_base, '.gz'];
input_info = niftiinfo(input_file);
input_info.PixelDimensions(1:3) = voxel_sizes_mm;
delete(input_file);
niftiwrite(input_data, input_base, input_info, 'Compressed', true);

mask_base = fullfile(case_dir, 'funcmask.nii');
niftiwrite(ones(volume_size, 'uint8'), mask_base, 'Compressed', true);
mask_file = [mask_base, '.gz'];

output_file = detrend_filter_smooth( ...
    input_file, mask_file, output_dir, fwhm_mm, [], 0);

if fwhm_mm == 0
    expected_output_file = fullfile(output_dir, 'synthetic_input.nii.gz');
    expected_volume = impulse;
else
    expected_output_file = fullfile(output_dir, ...
        [expected_tag, '_synthetic_input.nii.gz']);
    [expected_sigma_vox, ~, ~] = ...
        cicada_fwhm_mm_to_sigma_vox(fwhm_mm, voxel_sizes_mm);
    expected_volume = imgaussfilt3(impulse, expected_sigma_vox);
end

assert(strcmp(output_file, expected_output_file), ...
    'Unexpected output filename for case %s.', case_name);
assert(isfile(output_file), ...
    'Expected output file was not written for case %s.', case_name);

actual_data = niftiread(output_file);
expected_data = repmat(expected_volume, [1 1 1 n_timepoints]);
max_output_difference = max(abs( ...
    double(actual_data(:)) - double(expected_data(:))));
assert(max_output_difference < 1e-7, ...
    'Output values differed from direct expected smoothing for case %s.', ...
    case_name);

output_info = niftiinfo(output_file);
max_voxel_size_difference = max(abs( ...
    double(output_info.PixelDimensions(1:3)) - double(voxel_sizes_mm)));
assert(max_voxel_size_difference < 1e-6, ...
    'Output voxel metadata changed unexpectedly for case %s.', case_name);

if compare_legacy
    legacy_sigma_vox = ...
        (fwhm_mm / round(mean(voxel_sizes_mm))) / 2.355;
    legacy_volume = imgaussfilt3(impulse, legacy_sigma_vox);
    max_legacy_difference = max(abs( ...
        double(expected_volume(:)) - double(legacy_volume(:))));
    assert(max_legacy_difference < 2e-5, ...
        'Ordinary 2-mm output was not closely legacy-equivalent.');
end
end
