function CICADA_3_QC(cleaned_file, compare_file)
% Will compare QC of before and after CICADA, and compare to a compare file
% (8 parameter is standard). This will work for any cleaned_file, does not
% have to be cicada.
% cleaned dir is a directory within the output_dir given in the first two
% scripts. This will be either the cleaned_auto dir or the cleaned_manual
% dir
% While default compare file is 8p, you could feed it the cicada auto file
% if you were doing manual cicada step, for example. The file does need to
% already exist.

fprintf('\n')
close all

cleaned_file_info = dir(cleaned_file);
cleaned_dir = cleaned_file_info.folder;

if ~isfile(cleaned_file)
    fprintf(['Cannot find cleaned file at ', cleaned_file, '\n'])
    return;
end

cd(cleaned_dir)
cd('../')
[~, task_dir_name, ~]=fileparts(pwd); % grab task dir name
cd('../')
[~, ses_dir_name, ~]=fileparts(pwd); % grab ses dir name
cd('../')
[~, subj_dir_name, ~]=fileparts(pwd); % grab subj dir name
cd(cleaned_dir)


% Grab files, check if they exist 
orig_file_info = dir('*orig*.nii.gz'); % needs to be in cleaned_dir too
orig_file = [orig_file_info.folder, '/', orig_file_info.name];

% catch if there is no compare_file, and if so, do standard 8p (and auto
% later if cicada is manual version)
if ~exist('compare_file', 'var') || strcmp(compare_file, 'x')
    fprintf('Will compare to standard 8 parameter \n')
    compare_file_info = dir('*8p*');
    compare_file=[compare_file_info.folder, '/', compare_file_info.name]; % Give it default 8p to compare against
end


if ~isfile(compare_file)
    fprintf(['ERROR, cannot find compare file: ', compare_file, '\n'])
    return;
elseif ~isfile(orig_file)
    fprintf('ERROR, cannot find orig file: funcfile.nii.gz \n')
    return;
elseif ~isfile(cleaned_file)
    fprintf(['ERROR, cannot find Cleaned file: ', cleaned_file_info.name, '\n'])
    return;
end

orig_info = niftiinfo(orig_file);
tr = orig_info.PixelDimensions(4);

% check if comparing auto or manual, can base this off of cleaned file name
% then load relevant files from cicada
if contains(cleaned_file, '_CICADA_manual_')
    if isfile([cleaned_dir, '/../ic_manual_selection/DecisionVariables_Manual.mat'])
        load([cleaned_dir, '/../ic_manual_selection/DecisionVariables_Manual.mat']) %#ok<LOAD> 
        ic_select = 'manual';
        cicada_type = 'manual';
    else
        fprintf('Cannot find ../ic_manual_selection/DecisionVariables_Manual.mat \n')
        return
    end
elseif contains(cleaned_file, '_CICADA_auto_')
    if isfile([cleaned_dir, '/../ic_auto_selection/DecisionVariables_Auto.mat'])
        load([cleaned_dir, '/../ic_auto_selection/DecisionVariables_Auto.mat']) %#ok<LOAD> 
        ic_select = 'auto';
        cicada_type = 'auto';
    else
        fprintf('Cannot find ../ic_auto_selection/DecisionVariables_Auto.mat \n')
        return
    end
else
    % cleaned file is not CICADA
    cicada_type = 'none';
end


% calculate final estimation of DOF - probably want > 15 DOF, and/or at
% least 10% of percent_variance_kept from CICADA
% something to keep in mind is that you should actually calculate the power
% percentage being removed to be much more accurate - need to find the
% power fraction of bp from the data that has been cleaned by CICADA
% without frequency filtering thus far. Would use cicada_name.name without
% choice tag. Do later to improve this system.
if ~strcmp(cicada_type, 'none')
    DOF_estimate_final = Data.numvolumes .* Results.percent_variance_kept; %#ok<NODEF,USENS> 
    Results.num_ICs_kept = length(Results.signal_ICs);
    Results.num_ICs_total = length(Results.ICs);
    fprintf(['Estimate of Final DOF for ', cleaned_file_info.name, ' is %.2f\n'], DOF_estimate_final)
end



% Now, we can get the relevant correlation values and such!
fprintf(['Running comparison of ', cleaned_file_info.name, ' to ', dir(compare_file).name, '\n'])   

[denoised_Edge_GM_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_GM_corr, ...
    denoised_WMCSF_GM_corr, denoised_CSF_GM_corr, denoised_NotGM_GM_corr, denoised_GM_GM_autocorr, ...
    denoised_GM_mean] = CICADA_fileQC(cleaned_file, orig_file);

[compare_Edge_GM_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, compare_Outbrain_GM_corr, ...
    compare_WMCSF_GM_corr, compare_CSF_GM_corr, compare_NotGM_GM_corr, compare_GM_GM_autocorr, ...
    compare_GM_mean] = CICADA_fileQC(compare_file, orig_file);

[orig_Edge_GM_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, orig_Outbrain_GM_corr, ...
    orig_WMCSF_GM_corr, orig_CSF_GM_corr, orig_NotGM_GM_corr, orig_GM_GM_autocorr, ...
    orig_GM_mean] = CICADA_fileQC(orig_file, orig_file);


% Create QC folder (if it doesn't already exist)
if ~isfolder([cleaned_dir, '/../qc'])
    mkdir([cleaned_dir, '/../qc'])
end
cd([cleaned_dir, '/../qc'])

% extract naming tags
cleaned_tag = extractBetween(cleaned_file_info.name, [task_dir_name, '_'], '_bold.nii.gz');
cleaned_tag = cleaned_tag{:}; % need to expand it to just be char array
compare_tag = extractBetween(dir(compare_file).name, [task_dir_name, '_'], '_bold.nii.gz');
compare_tag = compare_tag{:}; % need to expand it to just be char array
orig_tag = 'orig';
prefix = [subj_dir_name, '_', ses_dir_name, '_task-', task_dir_name];

qc_naming = [prefix, '_', cleaned_tag, '_vs_', compare_tag]; 
qc_vals = [qc_naming, '_qc_vals.mat'];
qc_plots = [qc_naming, '_qc_plots.jpg'];

% save relevant plotting variables
fprintf('Saving Relevant QC Data\n')
title_string = [prefix, ': ', cleaned_tag, ', ', compare_tag, ', & ', orig_tag];
if ~strcmp(cicada_type, 'none')
    save(qc_vals, 'orig_NotGM_GM_corr', 'compare_NotGM_GM_corr', 'denoised_NotGM_GM_corr', ...
     'orig_GM_GM_autocorr', 'compare_GM_GM_autocorr', 'denoised_GM_GM_autocorr', ...
     'orig_Outbrain_GM_corr', 'compare_Outbrain_GM_corr', 'denoised_Outbrain_GM_corr', ...
     'orig_Edge_GM_corr', 'compare_Edge_GM_corr', 'denoised_Edge_GM_corr', ...
     'orig_WMCSF_GM_corr', 'compare_WMCSF_GM_corr', 'denoised_WMCSF_GM_corr', ...
     'orig_CSF_GM_corr', 'compare_CSF_GM_corr', 'denoised_CSF_GM_corr', ...
     'orig_DVARS_GM_corr', 'compare_DVARS_GM_corr', 'denoised_DVARS_GM_corr', ...
     'orig_FD_GM_corr', 'compare_FD_GM_corr', 'denoised_FD_GM_corr', ...
     'orig_GM_mean', 'compare_GM_mean', 'denoised_GM_mean', ...
     'title_string', 'orig_tag', 'compare_tag', 'cleaned_tag',...
     'cicada_type', 'ic_select', 'DOF_estimate_final', 'Results', 'tr', 'Data', 'Tables') %#ok<USENS> 
else
    save(qc_vals, 'orig_NotGM_GM_corr', 'compare_NotGM_GM_corr', 'denoised_NotGM_GM_corr', ...
     'orig_GM_GM_autocorr', 'compare_GM_GM_autocorr', 'denoised_GM_GM_autocorr', ...
     'orig_Outbrain_GM_corr', 'compare_Outbrain_GM_corr', 'denoised_Outbrain_GM_corr', ...
     'orig_Edge_GM_corr', 'compare_Edge_GM_corr', 'denoised_Edge_GM_corr', ...
     'orig_WMCSF_GM_corr', 'compare_WMCSF_GM_corr', 'denoised_WMCSF_GM_corr', ...
     'orig_CSF_GM_corr', 'compare_CSF_GM_corr', 'denoised_CSF_GM_corr', ...
     'orig_DVARS_GM_corr', 'compare_DVARS_GM_corr', 'denoised_DVARS_GM_corr', ...
     'orig_FD_GM_corr', 'compare_FD_GM_corr', 'denoised_FD_GM_corr', ...
     'orig_GM_mean', 'compare_GM_mean', 'denoised_GM_mean', ...
     'title_string', 'orig_tag', 'compare_tag', 'cleaned_tag',...
     'tr') %#ok<USENS> 
end


% Now, in the future, if you want to replot anything, you should have
% everything you need by loading the qc_values.mat file, and then run just
% the figure stuff below
denoised_name = 'CICADA';
plot_qc(denoised_Edge_GM_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_GM_corr, ...
    denoised_WMCSF_GM_corr, denoised_CSF_GM_corr, denoised_NotGM_GM_corr, denoised_GM_GM_autocorr, ...
    compare_Edge_GM_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, compare_Outbrain_GM_corr, ...
    compare_WMCSF_GM_corr, compare_CSF_GM_corr, compare_NotGM_GM_corr, compare_GM_GM_autocorr, ...
    orig_Edge_GM_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, orig_Outbrain_GM_corr, ...
    orig_WMCSF_GM_corr, orig_CSF_GM_corr, orig_NotGM_GM_corr, orig_GM_GM_autocorr, ...
    denoised_GM_mean, compare_GM_mean, orig_GM_mean, title_string, qc_plots, denoised_name)


fprintf('\n')
end