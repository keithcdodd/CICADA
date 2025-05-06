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

% if needed, get this function's filepath
current_path = mfilename('fullpath');
[qc_func_dir, ~, ~] = fileparts([current_path '.m']);
template_dir = [qc_func_dir, '/../templates'];

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

% get correlations, and make sure they are for the same voxel pairs for
% each method to allow for easy paired tests!
[denoised_Edge_Edge_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_Outbrain_corr, ...
    denoised_WMCSF_WMCSF_corr, denoised_CSF_CSF_corr, denoised_NotGM_NotGM_corr, denoised_GM_GM_corr, ...
    denoised_Suscept_Suscept_corr, denoised_GM_mean, denoised_mean_var_table, ...
    compare_Edge_Edge_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, compare_Outbrain_Outbrain_corr, ...
    compare_WMCSF_WMCSF_corr, compare_CSF_CSF_corr, compare_NotGM_NotGM_corr, compare_GM_GM_corr, ...
    compare_Suscept_Suscept_corr, compare_GM_mean, compare_mean_var_table, ...
    orig_Edge_Edge_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, orig_Outbrain_Outbrain_corr, ...
    orig_WMCSF_WMCSF_corr, orig_CSF_CSF_corr, orig_NotGM_NotGM_corr, orig_GM_GM_corr, ...
    orig_Suscept_Suscept_corr, orig_GM_mean, orig_mean_var_table] = CICADA_fileQC(cleaned_dir, cleaned_file, compare_file, orig_file);


% Calculate Denoising Success Parameter!
denoised_DS_prop = mean(abs(denoised_GM_GM_corr)) ./ mean(abs(denoised_NotGM_NotGM_corr));
compare_DS_prop = mean(abs(compare_GM_GM_corr)) ./ mean(abs(compare_NotGM_NotGM_corr));
orig_DS_prop = mean(abs(orig_GM_GM_corr)) ./ mean(abs(orig_NotGM_NotGM_corr));

denoised_DS = denoised_DS_prop / orig_DS_prop;
compare_DS = compare_DS_prop / orig_DS_prop;
orig_DS = 1; % always 1 since you compare to the original data


% Create QC folder (if it doesn't already exist)
if ~isfolder([cleaned_dir, '/../qc'])
    mkdir([cleaned_dir, '/../qc'])
end
cd([cleaned_dir, '/../qc'])
qc_output_dir = pwd;

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

% Paired nonparametric testing of mean magnitudes:
% denoise compare
[denoised_compare_Edge_Edge_p, ~, denoised_compare_Edge_Edge_stats] = signrank(abs(denoised_Edge_Edge_corr), abs(compare_Edge_Edge_corr));
[denoised_compare_FD_GM_p, ~, denoised_compare_FD_GM_stats] = signrank(abs(denoised_FD_GM_corr), abs(compare_FD_GM_corr));
[denoised_compare_DVARS_GM_p, ~, denoised_compare_DVARS_GM_stats] = signrank(abs(denoised_DVARS_GM_corr), abs(compare_DVARS_GM_corr));
[denoised_compare_Outbrain_Outbrain_p, ~, denoised_compare_Outbrain_Outbrain_stats] = signrank(abs(denoised_Outbrain_Outbrain_corr), abs(compare_Outbrain_Outbrain_corr));
[denoised_compare_WMCSF_WMCSF_p, ~, denoised_compare_WMCSF_WMCSF_stats] = signrank(abs(denoised_WMCSF_WMCSF_corr), abs(compare_WMCSF_WMCSF_corr));
[denoised_compare_CSF_CSF_p, ~, denoised_compare_CSF_CSF_stats] = signrank(abs(denoised_CSF_CSF_corr), abs(compare_CSF_CSF_corr));
[denoised_compare_Suscept_Suscept_p, ~, denoised_compare_Suscept_Suscept_stats] = signrank(abs(denoised_Suscept_Suscept_corr), abs(compare_Suscept_Suscept_corr));
[denoised_compare_NotGM_NotGM_p, ~, denoised_compare_NotGM_NotGM_stats] = signrank(abs(denoised_NotGM_NotGM_corr), abs(compare_NotGM_NotGM_corr));
[denoised_compare_GM_GM_p, ~, denoised_compare_GM_GM_stats] = signrank(abs(denoised_GM_GM_corr), abs(compare_GM_GM_corr));

% denoise orig
[denoised_orig_Edge_Edge_p, ~, denoised_orig_Edge_Edge_stats] = signrank(abs(denoised_Edge_Edge_corr), abs(orig_Edge_Edge_corr));
[denoised_orig_FD_GM_p, ~, denoised_orig_FD_GM_stats] = signrank(abs(denoised_FD_GM_corr), abs(orig_FD_GM_corr));
[denoised_orig_DVARS_GM_p, ~, denoised_orig_DVARS_GM_stats] = signrank(abs(denoised_DVARS_GM_corr), abs(orig_DVARS_GM_corr));
[denoised_orig_Outbrain_Outbrain_p, ~, denoised_orig_Outbrain_Outbrain_stats] = signrank(abs(denoised_Outbrain_Outbrain_corr), abs(orig_Outbrain_Outbrain_corr));
[denoised_orig_WMCSF_WMCSF_p, ~, denoised_orig_WMCSF_WMCSF_stats] = signrank(abs(denoised_WMCSF_WMCSF_corr), abs(orig_WMCSF_WMCSF_corr));
[denoised_orig_CSF_CSF_p, ~, denoised_orig_CSF_CSF_stats] = signrank(abs(denoised_CSF_CSF_corr), abs(orig_CSF_CSF_corr));
[denoised_orig_Suscept_Suscept_p, ~, denoised_orig_Suscept_Suscept_stats] = signrank(abs(denoised_Suscept_Suscept_corr), abs(orig_Suscept_Suscept_corr));
[denoised_orig_NotGM_NotGM_p, ~, denoised_orig_NotGM_NotGM_stats] = signrank(abs(denoised_NotGM_NotGM_corr), abs(orig_NotGM_NotGM_corr));
[denoised_orig_GM_GM_p, ~, denoised_orig_GM_GM_stats] = signrank(abs(denoised_GM_GM_corr), abs(orig_GM_GM_corr));


% combine into table
compare_significance_labels = {'cleaned_DS', 'compare_DS', 'cleaned_DS_compare_DS_ratio', ...
    'cleaned_Edge_mean_abs', 'compare_Edge_mean_abs', 'Edge_pval', ...
    'cleaned_FD_GM_mean_abs', 'compare_FD_GM_mean_abs', 'FD_GM_pval', ...
    'cleaned_DVARS_GM_mean_abs', 'compare_DVARS_GM_mean_abs', 'DVARS_GM_pval', ...
    'cleaned_Outbrain_mean_abs', 'compare_Outbrain_mean_abs', 'Outbrain_pval', ...
    'cleaned_WMCSF_mean_abs', 'compare_WMCSF_mean_abs', 'WMCSF_pval', ...
    'cleaned_CSF_mean_abs', 'compare_CSF_mean_abs', 'CSF_pval', ...
    'cleaned_Susceptibility_mean_abs', 'compare_Susceptibility_mean_abs', 'Susceptibility_pval', ...
    'cleaned_NotGM_mean_abs', 'compare_NotGM_mean_abs', 'NotGM_pval', ...
    'cleaned_GM_mean_abs', 'compare_GM_mean_abs', 'GM_pval'};

orig_significance_labels = {'cleaned_DS', 'orig_DS', 'cleaned_DS_orig_DS_ratio', ...
    'cleaned_Edge_mean_abs', 'orig_Edge_mean_abs', 'Edge_pval', ...
    'cleaned_FD_GM_mean_abs', 'orig_FD_GM_mean_abs', 'FD_GM_pval', ...
    'cleaned_DVARS_GM_mean_abs', 'orig_DVARS_GM_mean_abs', 'DVARS_GM_pval', ...
    'cleaned_Outbrain_mean_abs', 'orig_Outbrain_mean_abs', 'Outbrain_pval', ...
    'cleaned_WMCSF_mean_abs', 'orig_WMCSF_mean_abs', 'WMCSF_pval', ...
    'cleaned_CSF_mean_abs', 'orig_CSF_mean_abs', 'CSF_pval', ...
    'cleaned_Susceptibility_mean_abs', 'orig_Susceptibility_mean_abs', 'Susceptibility_pval', ...
    'cleaned_NotGM_mean_abs', 'orig_NotGM_mean_abs', 'NotGM_pval', ...
    'cleaned_GM_mean_abs', 'orig_GM_mean_abs', 'GM_pval'};

compare_significance_array = [denoised_DS, compare_DS, denoised_DS/compare_DS, mean(abs(denoised_Edge_Edge_corr)), mean(abs(compare_Edge_Edge_corr)), denoised_compare_Edge_Edge_p, ...
    mean(abs(denoised_FD_GM_corr)), mean(abs(compare_FD_GM_corr)), denoised_compare_FD_GM_p, ...
    mean(abs(denoised_DVARS_GM_corr)), mean(abs(compare_DVARS_GM_corr)), denoised_compare_DVARS_GM_p, ...
    mean(abs(denoised_Outbrain_Outbrain_corr)), mean(abs(compare_Outbrain_Outbrain_corr)), denoised_compare_Outbrain_Outbrain_p, ...
    mean(abs(denoised_WMCSF_WMCSF_corr)), mean(abs(compare_WMCSF_WMCSF_corr)), denoised_compare_WMCSF_WMCSF_p, ...
    mean(abs(denoised_CSF_CSF_corr)), mean(abs(compare_CSF_CSF_corr)), denoised_compare_CSF_CSF_p, ...
    mean(abs(denoised_Suscept_Suscept_corr)), mean(abs(compare_Suscept_Suscept_corr)), denoised_compare_Suscept_Suscept_p, ...
    mean(abs(denoised_NotGM_NotGM_corr)), mean(abs(compare_NotGM_NotGM_corr)), denoised_compare_NotGM_NotGM_p, ...
    mean(abs(denoised_GM_GM_corr)), mean(abs(compare_GM_GM_corr)), denoised_compare_GM_GM_p];

orig_significance_array = [denoised_DS, orig_DS, denoised_DS/orig_DS, mean(abs(denoised_Edge_Edge_corr)), mean(abs(orig_Edge_Edge_corr)), denoised_orig_Edge_Edge_p, ...
    mean(abs(denoised_FD_GM_corr)), mean(abs(orig_FD_GM_corr)), denoised_orig_FD_GM_p, ...
    mean(abs(denoised_DVARS_GM_corr)), mean(abs(orig_DVARS_GM_corr)), denoised_orig_DVARS_GM_p, ...
    mean(abs(denoised_Outbrain_Outbrain_corr)), mean(abs(orig_Outbrain_Outbrain_corr)), denoised_orig_Outbrain_Outbrain_p, ...
    mean(abs(denoised_WMCSF_WMCSF_corr)), mean(abs(orig_WMCSF_WMCSF_corr)), denoised_orig_WMCSF_WMCSF_p, ...
    mean(abs(denoised_CSF_CSF_corr)), mean(abs(orig_CSF_CSF_corr)), denoised_orig_CSF_CSF_p, ...
    mean(abs(denoised_Suscept_Suscept_corr)), mean(abs(orig_Suscept_Suscept_corr)), denoised_orig_Suscept_Suscept_p, ...
    mean(abs(denoised_NotGM_NotGM_corr)), mean(abs(orig_NotGM_NotGM_corr)), denoised_orig_NotGM_NotGM_p, ...
    mean(abs(denoised_GM_GM_corr)), mean(abs(orig_GM_GM_corr)), denoised_orig_GM_GM_p];

compare_significance_table = array2table(compare_significance_array, 'VariableNames', compare_significance_labels);
orig_significance_table = array2table(orig_significance_array, 'VariableNames', orig_significance_labels);

% write significance tables to qc folder: 
% signed rank of the mean of the absolute values of the sampled noise profile correlations 
% (paired, nonparametric, Wilcoxon)
writetable(compare_significance_table, [qc_output_dir, '/', prefix, '_', cleaned_tag, '_vs_', compare_tag, '_significance.csv'], 'Delimiter', ',')
writetable(orig_significance_table, [qc_output_dir, '/', prefix, '_', cleaned_tag, '_vs_', orig_tag, '_significance.csv'], 'Delimiter', ',')

delete tmp*.nii.gz

% Make Network Identifiability Image for QC (keep it within the best data with
% funcmask_constrained.nii.gz, remove outbrain because of potential high
% variance from noise)
%network_names = {"Visual", "Sensorimotor", "Dorsal Attention", "Salience", "Limbic_Reward", "Executive_Control", "Default"};
network_names = {"Default", "Sensorimotor", "Visual", "Salience", "Dorsal_Attention", "Executive_Control", "Limbic_Reward"};
network_identifiability_table = network_identifiability([template_dir, '/yeo_brainnetome_network_labels.nii.gz'], cleaned_file, compare_file, orig_file, [qc_output_dir, '/../funcmask.nii.gz'], [qc_output_dir, '/../region_masks/GM_mask.nii.gz'], 6, network_names, qc_output_dir);

% Create seed-based connectivity of the PCC to show how at least DMN
% connectivity improves!
%seed_based_conn(orig_file, [template_dir, '/PCC.nii.gz'], 3, [qc_output_dir, '/orig_PCC_conn'])
%seed_based_conn(compare_file, [template_dir, '/PCC.nii.gz'], 3, [qc_output_dir, '/compare_PCC_conn'])
%seed_based_conn(cleaned_file, [template_dir, '/PCC.nii.gz'], 3, [qc_output_dir, '/cicada_PCC_conn'])

% save relevant plotting variables
fprintf('Saving Relevant QC Data\n')
title_string = [prefix, ': ', cleaned_tag, ', ', compare_tag, ', & ', orig_tag];
if ~strcmp(cicada_type, 'none')
    save(qc_vals, 'orig_NotGM_NotGM_corr', 'compare_NotGM_NotGM_corr', 'denoised_NotGM_NotGM_corr', ...
     'orig_GM_GM_corr', 'compare_GM_GM_corr', 'denoised_GM_GM_corr', ...
     'orig_Outbrain_Outbrain_corr', 'compare_Outbrain_Outbrain_corr', 'denoised_Outbrain_Outbrain_corr', ...
     'orig_Edge_Edge_corr', 'compare_Edge_Edge_corr', 'denoised_Edge_Edge_corr', ...
     'orig_WMCSF_WMCSF_corr', 'compare_WMCSF_WMCSF_corr', 'denoised_WMCSF_WMCSF_corr', ...
     'orig_CSF_CSF_corr', 'compare_CSF_CSF_corr', 'denoised_CSF_CSF_corr', ...
     'orig_Suscept_Suscept_corr', 'compare_Suscept_Suscept_corr', 'denoised_Suscept_Suscept_corr', ...
     'orig_DVARS_GM_corr', 'compare_DVARS_GM_corr', 'denoised_DVARS_GM_corr', ...
     'orig_FD_GM_corr', 'compare_FD_GM_corr', 'denoised_FD_GM_corr', ...
     'compare_significance_table', 'orig_significance_table', ...
     'orig_DS_prop', 'compare_DS_prop', 'denoised_DS_prop', ...
     'orig_DS', 'compare_DS', 'denoised_DS', 'network_identifiability_table', ...
     'orig_GM_mean', 'compare_GM_mean', 'denoised_GM_mean', ...
     'orig_mean_var_table', 'compare_mean_var_table', 'denoised_mean_var_table', ...
     'title_string', 'orig_tag', 'compare_tag', 'cleaned_tag',...
     'cicada_type', 'ic_select', 'DOF_estimate_final', 'Results', 'tr', 'Data', 'Tables') %#ok<USENS> 
else
    save(qc_vals, 'orig_NotGM_NotGM_corr', 'compare_NotGM_NotGM_corr', 'denoised_NotGM_NotGM_corr', ...
     'orig_GM_GM_corr', 'compare_GM_GM_corr', 'denoised_GM_GM_corr', ...
     'orig_Outbrain_Outbrain_corr', 'compare_Outbrain_Outbrain_corr', 'denoised_Outbrain_Outbrain_corr', ...
     'orig_Edge_Edge_corr', 'compare_Edge_Edge_corr', 'denoised_Edge_Edge_corr', ...
     'orig_WMCSF_WMCSF_corr', 'compare_WMCSF_WMCSF_corr', 'denoised_WMCSF_WMCSF_corr', ...
     'orig_CSF_CSF_corr', 'compare_CSF_CSF_corr', 'denoised_CSF_CSF_corr', ...
     'orig_Suscept_Suscept_corr', 'compare_Suscept_Suscept_corr', 'denoised_Suscept_Suscept_corr', ...
     'orig_DVARS_GM_corr', 'compare_DVARS_GM_corr', 'denoised_DVARS_GM_corr', ...
     'orig_FD_GM_corr', 'compare_FD_GM_corr', 'denoised_FD_GM_corr', ...
     'compare_significance_table', 'orig_significance_table', ...
     'orig_DS_prop', 'compare_DS_prop', 'denoised_DS_prop', ...
     'orig_DS', 'compare_DS', 'denoised_DS', 'network_identifiability_table', ...
     'orig_GM_mean', 'compare_GM_mean', 'denoised_GM_mean', ...
     'orig_mean_var_table', 'compare_mean_var_table', 'denoised_mean_var_table', ...
     'title_string', 'orig_tag', 'compare_tag', 'cleaned_tag',...
     'tr')  
end




% Now, in the future, if you want to replot anything, you should have
% everything you need by loading the qc_values.mat file, and then run just
% the figure stuff below
fprintf('Plotting QC\n')
denoised_name = 'CICADA';
if strcmp(cicada_type, 'none')
    % if it is not CICADA, label the cleaned file more appropriately
    denoised_name = cleaned_tag;
end
plot_qc(denoised_Edge_Edge_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_Outbrain_corr, ...
    denoised_WMCSF_WMCSF_corr, denoised_CSF_CSF_corr, denoised_NotGM_NotGM_corr, denoised_GM_GM_corr, denoised_Suscept_Suscept_corr, ...
    compare_Edge_Edge_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, compare_Outbrain_Outbrain_corr, ...
    compare_WMCSF_WMCSF_corr, compare_CSF_CSF_corr, compare_NotGM_NotGM_corr, compare_GM_GM_corr, compare_Suscept_Suscept_corr, ...
    orig_Edge_Edge_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, orig_Outbrain_Outbrain_corr, ...
    orig_WMCSF_WMCSF_corr, orig_CSF_CSF_corr, orig_NotGM_NotGM_corr, orig_GM_GM_corr, orig_Suscept_Suscept_corr, ...
    denoised_GM_mean, compare_GM_mean, orig_GM_mean, title_string, qc_plots, denoised_name)


fprintf('\n')
end