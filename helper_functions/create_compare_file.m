function compare_file = create_compare_file(output_dir, compare_file)
% will look for, and try to create, compare_file if it does not exist
% valid tags like {'8p', '9p'} etc.
% point is to look for a compare file given the compare_tag, and if not
% there, try to make it. This is only relevant if basescript 2 has already
% been run! Relies on finding orig data!

cleaned_dir = [output_dir, '/cleaned'];
task_dir = output_dir;

compare_file_info = dir([cleaned_dir, '/*', compare_file, '*']); % help see if file exists as it should

if isempty(compare_file_info)
    % file does not exist yet, try to make it!
    fprintf('Cannot find the compare file in cleaned directory... will try to recreate it from files from basescript 2!\n')
    compare_tag = compare_file; % since it is truly a tag here, not a full filepath

    % make sure the regressors file exists!
    if isfile([task_dir, '/regressors_timeseries/', compare_tag, '_regressors_intercept.mat'])

        % great, we should be able to make it!
        % can get the appropriate prefix and suffix from the orig data (which SHOULD exist):
        orig_file_info = dir([cleaned_dir, '/*', compare_file, '*']);
        parts = split(orig_file_info.name, '_orig_');
        prefix = parts{1};
        suffix = parts{2};

        % grab Tmean to add back in
        tmean_command = ['fslmaths ', task_dir, '/funcfile.nii.gz -Tmean ', task_dir, '/tmean_funcfile.nii.gz'];
        [~, ~] = call_fsl(tmean_command);
        
        % Now, run regression for comparison file of interest!
        compare_regress_command = ['fsl_glm -i ', task_dir, '/funcfile.nii.gz -d ', task_dir, ...
            '/regressors_timeseries/', compare_tag, '_regressors_intercept.mat -m ', task_dir, '/funcmask.nii.gz ', ...
            '--out_res=', cleaned_dir, '/', prefix, '_', compare_tag, '_', suffix];
        fprintf(['Running: ', compare_regress_command, '\n'])
        [~, ~] = call_fsl(compare_regress_command);
        
        % Because fsl_glm is dumb, it resets the TR to 1... need to run fslmerge
        % with tr option to reset the tr correctly...
        reset_tr_compare_command = ['fslmerge -tr ', cleaned_dir, '/', prefix, '_', compare_tag, '_', suffix, ' ', cleaned_dir, '/', prefix, '_', compare_tag, '_', suffix, ' ', num2str(TR)];
        [~, ~] = call_fsl(reset_tr_compare_command);
        
        % add Tmean back onto compare file 
        tmean_add_compare_command = ['fslmaths ', cleaned_dir, '/', prefix, '_', compare_tag, '_', suffix, ' -add ', task_dir, '/tmean_funcfile.nii.gz ', cleaned_dir, '/', prefix, '_', compare_tag, '_', suffix];
        [~, ~] = call_fsl(tmean_add_compare_command);

        compare_file_info = dir([cleaned_dir, '/*', compare_tag, '*']); % Nowwwww, it should find it!
    else
        fprintf('Missing (?) regressor mat file for fsl_glm that should have been made by the 2nd basescript... that is bad...\n')
        fprintf('Will just use 8p which should always be made by default!\n')
        compare_file_info = dir([cleaned_dir, '/*8p*']); % 8 parameter! Should always exist since it is run by default!
    end

    
end

compare_file = [compare_file_info.folder, '/', compare_file_info.name]; % Update it to match to the tagged type of comparison denoiser

if ~isfile(compare_file)
    fprintf('Something is probably messed up... could not find, create, or even use default 8p compare file...?? \nA pathing issue maybe, or you tried to run this before Auto CICADA?\n')
    compare_file = '';
end

end

        


