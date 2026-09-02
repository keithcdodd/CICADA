function test_output_smoothing_geometry
%TEST_OUTPUT_SMOOTHING_GEOMETRY Numeric tests for physical output smoothing.

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));
addpath(fullfile(repo_root, 'helper_functions'));

fwhm_mm = 6;
exact_factor = 2 * sqrt(2 * log(2));

% Ordinary 2-mm isotropic data remains numerically equivalent to the
% historical scalar conversion (which used the rounded constant 2.355).
[sigma_iso, sigma_mm, vox_iso] = ...
    cicada_fwhm_mm_to_sigma_vox(fwhm_mm, [2 2 2]);
legacy_sigma = (fwhm_mm / 2) / 2.355;
assert(all(abs(sigma_iso - legacy_sigma) < 1e-4));
assert(abs(sigma_mm - fwhm_mm / exact_factor) < 1e-12);
assert(isequal(vox_iso, [2 2 2]));

% Per-axis conversion preserves one physical FWHM for anisotropic voxels.
[sigma_aniso, ~, vox_aniso] = ...
    cicada_fwhm_mm_to_sigma_vox(fwhm_mm, [2 2 3]);
recovered_fwhm = sigma_aniso .* vox_aniso .* exact_factor;
assert(all(abs(recovered_fwhm - fwhm_mm) < 1e-12));

% Non-integer voxel dimensions are not rounded.
voxel_sizes = [2.4 2.4 2.7];
[sigma_noninteger, ~, returned_voxels] = ...
    cicada_fwhm_mm_to_sigma_vox(5.5, voxel_sizes);
assert(isequal(returned_voxels, voxel_sizes));
assert(all(abs(sigma_noninteger .* voxel_sizes .* exact_factor - 5.5) < 1e-12));

% Integer names remain compatible; non-integer names remain distinct.
assert(strcmp(cicada_smoothing_fwhm_tag(6), 's6'));
assert(strcmp(cicada_smoothing_fwhm_tag(6.4), 's6p4'));
assert(~strcmp(cicada_smoothing_fwhm_tag(6.4), ...
    cicada_smoothing_fwhm_tag(6.49)));

assert_throws(@() cicada_fwhm_mm_to_sigma_vox(-1, [2 2 2]));
assert_throws(@() cicada_fwhm_mm_to_sigma_vox(NaN, [2 2 2]));
assert_throws(@() cicada_fwhm_mm_to_sigma_vox(6, [2 0 2]));
assert_throws(@() cicada_fwhm_mm_to_sigma_vox(6, [2 Inf 2]));
assert_throws(@() cicada_fwhm_mm_to_sigma_vox(6, [2 2]));

fprintf('PASS: physical output-smoothing geometry tests.\n');
end

function assert_throws(test_function)
did_throw = false;
try
    test_function();
catch
    did_throw = true;
end
assert(did_throw, 'Expected validation to reject invalid input.');
end
