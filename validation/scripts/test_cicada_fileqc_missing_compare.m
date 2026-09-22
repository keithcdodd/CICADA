%% Regression test: CICADA_fileQC catches a missing comparison file
%
% The denoised file exists, while compare_file deliberately does not.
% CICADA_fileQC should stop at the comparison-file validation check and
% report the correct missing path.

script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(script_dir));

addpath(fullfile(repo_root, 'helper_functions'));

test_dir = tempname;
mkdir(test_dir);
cleanup_obj = onCleanup(@() rmdir(test_dir, 's'));

denoised_file = fullfile(test_dir, 'denoised_dummy.nii.gz');
compare_file  = fullfile(test_dir, 'missing_compare.nii.gz');
orig_file     = fullfile(test_dir, 'orig_dummy.nii.gz');

% CICADA_fileQC only needs denoised_file to exist before reaching
% the compare-file check, so an empty placeholder is sufficient.
fid = fopen(denoised_file, 'w');
assert(fid ~= -1, 'Could not create temporary denoised-file placeholder.');
fclose(fid);

assert(isfile(denoised_file));
assert(~isfile(compare_file));

captured = evalc( ...
    'CICADA_fileQC(test_dir, denoised_file, compare_file, orig_file);');

assert(contains(captured, 'Cannot find compare file at'), ...
    'CICADA_fileQC did not report a missing comparison file.');

assert(contains(captured, compare_file), ...
    'CICADA_fileQC reported the wrong comparison-file path.');

assert(~contains(captured, 'Cannot find denoised file'), ...
    'CICADA_fileQC incorrectly failed the denoised-file check.');

fprintf('PASS: CICADA_fileQC correctly detects a missing comparison file.\n');