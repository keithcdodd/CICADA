function compare_file = create_compare_file(output_dir, compare_file)
% CREATE_COMPARE_FILE Recreate a missing standard comparison denoising file.
%
% This is primarily used after CICADA_2_AutoLabeling has already generated
% the standard confound-regression design matrices.
%
% compare_file may be:
%   1. An existing comparison filepath -> returned unchanged, or
%   2. A comparison tag such as '12p' -> recreated if possible.
%
% If the requested regression design is unavailable, CICADA falls back to
% the standard 8p comparison when possible.

task_dir = output_dir;
cleaned_dir = fullfile(task_dir, 'cleaned');

if isfile(compare_file)
    return
end

if nargin < 2 || isempty(compare_file) || ~ischar(compare_file)
    compare_tag = '8p';
else
    compare_tag = compare_file;
end


%% Check whether the requested comparison already exists

compare_file_info = ...
    dir(fullfile(cleaned_dir, ['*_', compare_tag, '_*.nii.gz']));

if isscalar(compare_file_info)

    compare_file = fullfile( ...
        compare_file_info(1).folder, ...
        compare_file_info(1).name);

    return

elseif numel(compare_file_info) > 1

    error('CICADA:AmbiguousCompareFile', ...
        'Multiple cleaned comparison files were found for tag "%s" in:\n%s', ...
        compare_tag, cleaned_dir);

end


%% Find the required regression design

regressor_file = fullfile( ...
    task_dir, ...
    'regressors_timeseries', ...
    [compare_tag, '_regressors_intercept.mat']);

if ~isfile(regressor_file)

    if ~strcmp(compare_tag, '8p')

        warning('CICADA:MissingComparisonRegressors', ...
            ['Cannot find regressors for requested comparison "%s". ' ...
             'Falling back to standard 8p.'], ...
            compare_tag);

        compare_tag = '8p';

        % First see whether 8p already exists.
        compare_file_info = ...
            dir(fullfile(cleaned_dir, '*_8p_*.nii.gz'));

        if isscalar(compare_file_info)

            compare_file = fullfile( ...
                compare_file_info(1).folder, ...
                compare_file_info(1).name);

            return

        elseif numel(compare_file_info) > 1

            error('CICADA:AmbiguousCompareFile', ...
                'Multiple cleaned 8p comparison files were found in:\n%s', ...
                cleaned_dir);

        end

        regressor_file = fullfile( ...
            task_dir, ...
            'regressors_timeseries', ...
            '8p_regressors_intercept.mat');

    end

end

if ~isfile(regressor_file)

    warning('CICADA:ComparisonRegressorsUnavailable', ...
        ['Cannot recreate comparison "%s" because the required regression ' ...
         'design matrix is missing:\n%s'], ...
        compare_tag, regressor_file);

    compare_file = '';
    return

end


%% Find the original functional-data copy

orig_file_info = ...
    dir(fullfile(cleaned_dir, '*_orig_*.nii.gz'));

if isempty(orig_file_info)

    warning('CICADA:OriginalComparisonTemplateMissing', ...
        ['Cannot recreate comparison "%s" because no *_orig_*.nii.gz ' ...
         'file exists in:\n%s'], ...
        compare_tag, cleaned_dir);

    compare_file = '';
    return

elseif numel(orig_file_info) > 1

    error('CICADA:AmbiguousOriginalFile', ...
        ['Multiple *_orig_*.nii.gz files were found in:\n%s\n' ...
         'Cannot determine the comparison filename safely.'], ...
        cleaned_dir);

end


%% Required functional inputs

funcfile = fullfile(task_dir, 'funcfile.nii.gz');
funcmask = fullfile(task_dir, 'funcmask.nii.gz');

if ~isfile(funcfile)

    warning('CICADA:MissingFuncfile', ...
        'Cannot recreate comparison because funcfile.nii.gz is missing.');

    compare_file = '';
    return

end

if ~isfile(funcmask)

    warning('CICADA:MissingFuncmask', ...
        'Cannot recreate comparison because funcmask.nii.gz is missing.');

    compare_file = '';
    return

end


%% Derive comparison filename from the existing original filename

orig_name = orig_file_info(1).name;

compare_name = strrep( ...
    orig_name, ...
    '_orig_', ...
    ['_', compare_tag, '_']);

compare_file = fullfile(cleaned_dir, compare_name);


%% Read the actual TR from the functional image

func_info = niftiinfo(funcfile);

if numel(func_info.PixelDimensions) < 4

    error('CICADA:MissingFunctionalTR', ...
        'Could not determine the functional TR from funcfile.nii.gz.');

end

TR = func_info.PixelDimensions(4);


%% Calculate temporal mean

tmean_file = fullfile(task_dir, 'tmean_funcfile.nii.gz');

tmean_command = sprintf( ...
    'fslmaths "%s" -Tmean "%s"', ...
    funcfile, tmean_file);

[status, msg] = call_fsl(tmean_command);

if status ~= 0
    error('CICADA:FSLFailure', ...
        'Failed to calculate functional temporal mean: %s', string(msg));
end


%% Recreate comparison residuals

fprintf('Recreating %s comparison file.\n', compare_tag)

compare_regress_command = sprintf( ...
    'fsl_glm -i "%s" -d "%s" -m "%s" --out_res="%s"', ...
    funcfile, regressor_file, funcmask, compare_file);

fprintf('Running: %s\n', compare_regress_command)

[status, msg] = call_fsl(compare_regress_command);

if status ~= 0
    error('CICADA:FSLFailure', ...
        'Failed to recreate %s comparison: %s', ...
        compare_tag, string(msg));
end


%% Restore TR

% fsl_glm residual output resets TR to 1 second, so restore the original TR.
reset_tr_command = sprintf( ...
    'fslmerge -tr "%s" "%s" %.15g', ...
    compare_file, compare_file, TR);

[status, msg] = call_fsl(reset_tr_command);

if status ~= 0
    error('CICADA:FSLFailure', ...
        'Failed to restore TR for %s comparison: %s', ...
        compare_tag, string(msg));
end


%% Add temporal mean back

tmean_add_command = sprintf( ...
    'fslmaths "%s" -add "%s" "%s"', ...
    compare_file, tmean_file, compare_file);

[status, msg] = call_fsl(tmean_add_command);

if status ~= 0
    error('CICADA:FSLFailure', ...
        'Failed to restore temporal mean for %s comparison: %s', ...
        compare_tag, string(msg));
end


%% Final validation

if ~isfile(compare_file)

    warning('CICADA:ComparisonCreationFailed', ...
        ['CICADA attempted to create comparison "%s", but the expected ' ...
         'output file was not found:\n%s'], ...
        compare_tag, compare_file);

    compare_file = '';

    return

end

fprintf('Created comparison file:\n  %s\n', compare_file)

end