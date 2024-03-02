function CICADA_3_QC(cleaned_dir, compare_file)
% Will compare QC of before and after CICADA, and compare to a compare file
% (8 parameter is standard)
% cleaned dir is a directory within the output_dir given in the first two
% scripts. This will be either the cleaned_auto dir or the cleaned_manual
% dir
% While default compare file is 8p, you could feed it the cicada auto file
% if you were doing manual cicada step, for example. The file does need to
% already exist.

fprintf('\n')
close all

cd(cleaned_dir)
cd('../')
[~, task_dir_name, ~]=fileparts(pwd); % grab task dir name
cd('../')
[~, ses_dir_name, ~]=fileparts(pwd); % grab ses dir name
cd('../')
[~, subj_dir_name, ~]=fileparts(pwd); % grab subj dir name
cd(cleaned_dir)

% need to check if manual or auto, because there will be two cicada files
% if manual
[~, cleaned_dir_name, ~]=fileparts(pwd);

if contains(cleaned_dir_name, 'manual')
    cicada_type = 'manual';
else
    cicada_type = 'auto';
end

% Grab files, check if they exist 
cicada_file_info = dir(['*CICADA*', cicada_type, '*nonagg*.nii.gz']); % grabs the auto if cleaned auto, grabs manual if cleaned manual
cicada_file = [cicada_file_info.folder, '/', cicada_file_info.name];
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
elseif ~isfile(cicada_file)
    fprintf(['ERROR, cannot find CICADA file: ', cicada_file_info.name, '\n'])
    return;
end

orig_info = niftiinfo(orig_file);
tr = orig_info.PixelDimensions(4);

% check if comparing auto or manual, can base this off of cleaned dir name
% then load relevant files from cicada
if contains(cleaned_dir_name, 'manual')
    if isfile([cleaned_dir, '/../ic_manual_selection/DecisionVariables_Manual.mat'])
        load([cleaned_dir, '/../ic_manual_selection/DecisionVariables_Manual.mat']) %#ok<LOAD> 
        ic_select = 'manual';
    else
        fprintf('Cannot find ../ic_manual_selection/DecisionVariables_Manual.mat \n')
        return
    end
else
    if isfile([cleaned_dir, '/../ic_auto_selection/DecisionVariables_Auto.mat'])
        load([cleaned_dir, '/../ic_auto_selection/DecisionVariables_Auto.mat']) %#ok<LOAD> 
        ic_select = 'auto';
    else
        fprintf('Cannot find ../ic_auto_selection/DecisionVariables_Auto.mat \n')
        return
    end
end


% calculate final estimation of DOF - probably want > 15 DOF, and/or at
% least 10% of percent_variance_kept from CICADA
% something to keep in mind is that you should actually calculate the power
% percentage being removed to be much more accurate - need to find the
% power fraction of bp from the data that has been cleaned by CICADA
% without frequency filtering thus far. Would use cicada_name.name without
% choice tag. Do later to improve this system.
DOF_estimate_final = Data.numvolumes .* Results.percent_variance_kept; %#ok<NODEF,USENS> 
Results.num_ICs_kept = length(Results.signal_ICs);
Results.num_ICs_total = length(Results.ICs);
fprintf(['Estimate of Final DOF for ', cicada_file_info.name, ' is %.2f\n'], DOF_estimate_final)


% Now, we can get the relevant correlation values and such!
fprintf(['Running comparison of ', cicada_file_info.name, ' to ', dir(compare_file).name, '\n'])   

[denoised_Edge_GM_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_GM_corr, ...
    denoised_WMCSF_GM_corr, denoised_CSF_GM_corr, denoised_NotGM_GM_corr, denoised_GM_GM_autocorr, ...
    denoised_GM_mean] = CICADA_fileQC(cicada_file, orig_file);

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
delete([cleaned_dir, '/../qc/*']) % delete old files to save potential space
cd([cleaned_dir, '/../qc'])

% extract naming tags
cicada_tag = extractBetween(cicada_file_info.name, [task_dir_name, '_'], '_bold.nii.gz');
cicada_tag = cicada_tag{:}; % need to expand it to just be char array
compare_tag = extractBetween(dir(compare_file).name, [task_dir_name, '_'], '_bold.nii.gz');
compare_tag = compare_tag{:}; % need to expand it to just be char array
orig_tag = 'orig';
prefix = [subj_dir_name, '_', ses_dir_name, '_task-', task_dir_name];

qc_naming = [prefix, '_', cicada_tag, '_vs_', compare_tag]; 
qc_vals = [qc_naming, '_qc_vals.mat'];
qc_plots = [qc_naming, '_qc_plots.jpg'];

% save relevant plotting variables
fprintf('Saving Relevant QC Data\n')
title_string = [prefix, ': ', cicada_tag, ', ', compare_tag, ', & ', orig_tag];
save('qc_vals', 'orig_NotGM_GM_corr', 'compare_NotGM_GM_corr', 'denoised_NotGM_GM_corr', ...
     'orig_GM_GM_autocorr', 'compare_GM_GM_autocorr', 'denoised_GM_GM_autocorr', ...
     'orig_Outbrain_GM_corr', 'compare_Outbrain_GM_corr', 'denoised_Outbrain_GM_corr', ...
     'orig_Edge_GM_corr', 'compare_Edge_GM_corr', 'denoised_Edge_GM_corr', ...
     'orig_WMCSF_GM_corr', 'compare_WMCSF_GM_corr', 'denoised_WMCSF_GM_corr', ...
     'orig_CSF_GM_corr', 'compare_CSF_GM_corr', 'denoised_CSF_GM_corr', ...
     'orig_DVARS_GM_corr', 'compare_DVARS_GM_corr', 'denoised_DVARS_GM_corr', ...
     'orig_FD_GM_corr', 'compare_FD_GM_corr', 'denoised_FD_GM_corr', ...
     'denoised_GM_mean', 'compare_GM_mean', 'orig_GM_mean', ...
     'title_string', 'orig_tag', 'compare_tag', 'cicada_tag',...
     'ic_select', 'DOF_estimate_final', 'Results', 'tr', 'Data', 'Tables') %#ok<USENS> 

% Now, in the future, if you want to replot anything, you should have
% everything you need by loading the qc_values.mat file, and then run just
% the figure stuff below

fprintf('Creating Figures\n')
figure('Position', [50, 50, 2000, 900])
t = tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');
title(t, title_string, 'Interpreter', 'none')

% GM_Edge corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('Edge-GM Corr ', 'Interpreter', 'none')
histogram(orig_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_fd corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against overall movement.
nexttile
hold on
title('FD-GM Corr ', 'Interpreter', 'none')
histogram(orig_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_dvars corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against global signal changes (including large
% movement)
nexttile
hold on
title('DVARS-GM Corr ', 'Interpreter', 'none')
histogram(orig_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_Outbrain corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('Outbrain-GM Corr ', 'Interpreter', 'none')
histogram(orig_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_WMCSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('WMCSF-GM Corr ', 'Interpreter', 'none')
histogram(orig_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_CSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('CSF-GM Corr ', 'Interpreter', 'none')
histogram(orig_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_NotGM corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that. Want to see the denoised and orig be significantly different
nexttile
hold on
title('NotGM-GM Corr ', 'Interpreter', 'none')
histogram(orig_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
hold off

% GM_GM autocorr to orig: Want to see signal correlations maintained, in comparison to
% GM_NotGM corr, but still reduced (removed erroneous/global correlations).
% get color order right to match other plots
colorord = get(gca, 'colororder');

nexttile
hold on
title('GM-GM Auto Corr to Orig', 'Interpreter', 'none')
histogram(compare_GM_GM_autocorr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
histogram(denoised_GM_GM_autocorr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
legend('compare', 'CICADA', 'Interpreter', 'none')
hold off

% Plot GS changes before and after Denoising - it is easier to
% detrend it all first for comparison
nexttile
hold on
title('Mean GM Signal (Detrended)', 'Interpreter', 'none')
plot(detrend(orig_GM_mean,2), 'LineWidth', 1.5)
plot(detrend(compare_GM_mean,2), 'LineWidth', 1.5)
plot(detrend(denoised_GM_mean,2), 'LineWidth', 1.5)
legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
xlim([0, length(orig_GM_mean)])
hold off


fprintf('Saving Figure\n')
exportgraphics(t, qc_plots, 'Resolution', 300)

fprintf('\n')
end