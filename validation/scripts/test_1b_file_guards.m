function test_1b_file_guards
% Public-safe functional checks for early CICADA_3_QC file guards.

repo_root = fileparts(fileparts(fileparts(mfilename('fullpath'))));
addpath(fullfile(repo_root, 'basescripts'));
addpath(fullfile(repo_root, 'helper_functions'));

original_dir = pwd;
test_root = tempname;
mkdir(test_root);
cleanup_obj = onCleanup(@() cleanup_test_root(test_root, original_dir)); %#ok<NASGU>

cleaned_dir = fullfile(test_root, 'sub-test', 'ses-test', 'task-test', 'cleaned');
mkdir(cleaned_dir);

% 1. Missing cleaned file must return before dir() dereference.
missing_cleaned = fullfile(cleaned_dir, 'missing_cleaned.nii.gz');
out = evalc("CICADA_3_QC(missing_cleaned, '8p');");
assert(contains(out, 'Cannot find cleaned file'));

% Placeholder cleaned file. Later cases must still return before NIfTI read.
cleaned_file = fullfile(cleaned_dir, ...
    'sub-test_ses-test_task-test_CICADA_auto_bold.nii.gz');
touch_file(cleaned_file);

% 2. Missing original file must fail cleanly.
out = evalc("CICADA_3_QC(cleaned_file, '8p');");
assert(contains(out, 'cannot find orig file matching'));

% 3. Multiple original files must be rejected explicitly.
orig1 = fullfile(cleaned_dir, 'sub-test_task-test_orig_bold.nii.gz');
orig2 = fullfile(cleaned_dir, 'sub-test_task-test_orig_copy_bold.nii.gz');
touch_file(orig1);
touch_file(orig2);

out = evalc("CICADA_3_QC(cleaned_file, '8p');");
assert(contains(out, 'found multiple orig files'));

delete(orig2);

% 4. One original file but no comparison file must fail cleanly
% before any NIfTI read.
out = evalc("CICADA_3_QC(cleaned_file, '8p');");
assert(contains(out, 'cannot find compare file'));

fprintf('PASS: CICADA_3_QC early file guards.\n');
end

function touch_file(file_path)
fid = fopen(file_path, 'w');
assert(fid ~= -1, 'Could not create temporary test file.');
fclose(fid);
end

function cleanup_test_root(test_root, original_dir)
% CICADA_3_QC changes directories internally, so restore a persistent
% directory before removing the temporary test tree.
if isfolder(original_dir)
    cd(original_dir);
end
if isfolder(test_root)
    rmdir(test_root, 's');
end
end
