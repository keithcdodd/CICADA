function CADICA_4_QC(cadicafol, participant_id, sess_id, task_id, prefix, suffix, compare_file_tag, ic_select, aggressionchoice, filterchoice, smoothchoice, lowfreq_cutoff, highfreq_cutoff, repetitiontime)
% CADICA_4_QC: To be run after CADICA Clean
% compare_filename: the "raw data" you want to compare to the denoised data.
    % Example options are 8p_funcfile_bold.nii.gz, or
    % 9p_funcfile_bold.nii.gz. But this can be altered as one sees fit in
    % this function. Can input smoothed and/or filtered versions. If it is
    % filtered, it will also be smoothed. A good default is to compare
    % files before smoothing or filtering, so one can focus on the effects
    % of the different regression styles.
    % ic_select: manual or auto, default is manual
% aggressionchoice: agg or nonagg (which type of ICA denoising are you
    % comparing to?). Default is nonagg (optional)
% filterchoice: bp, hp, or nf (which type of filtering are you comparingf
    % to, bandpass, highpass, or no filter?). Default is nf (optional)
% smoothchoice: s or ns (smoothing or no smoothing) for comparison. Either
% option is good. To compare just the raw effects of regression, one could
% pick ns (no smoothing). Default is no smoothing. If filter choice is hp
% or bd (not no filter), then smoothing will be added on.
% lowfreq_cutoff: In Hz. If not given, then will be set to 0 if
% filterchoice is nf, and 0.008 if filterchoice is hp or bp. (optional)
% highfreq_cutoff: In Hz. If not given, then will be set to 1/(2*tr) or
% 0.15 if filterchoice is bp
% repetitiontime: tr (needs to calculate estimated Degrees of Freedom)

% e.g. CADICA_4QC('/path/cadica_MNI152NLin6Asym_res-02', '102', '01',
% 'foodrun_run-01', 'funcfile_bold.nii.gz', 'nonagg', 'bp', 0.008, 0.15)

% the most important things to see include reduced DVARS and FD correlation
% (movement is reduced), and global signal should show less severe spiking 
% (unless the comparison denoising was overly aggressive e.g., 9p potentially).
% The correlation between GM and nonGM should also be reduced a little.

% Note: All comparisons are to the original raw data's image

fprintf('\n')
close all

tr=repetitiontime;

% get to relevant folder
cd(cadicafol)
subj_fol = ['sub-', participant_id];
cd(subj_fol)
sess_fol = ['ses-', sess_id];
cd(sess_fol)
cd(task_id)

% check optional inputs
if ~exist('filterchoice', 'var')
    filterchoice='nf';
end

if ~exist('smoothchoice', 'var')
    smoothchoice='ns';
end

if ~exist('aggressionchoice', 'var')
    aggressionchoice='nonagg';
end

if ~strcmp(aggressionchoice, 'agg')
    aggressionchoice='nonagg';
end

if ~strcmp(ic_select, 'auto')
    ic_select = 'manual';
    if isfile('./ic_manual_selection/DecisionVariables_Manual.mat')
        load('./ic_manual_selection/DecisionVariables_Manual.mat')
    else
        fprintf('Cannot find ./ic_manual_selection/DecisionVariables_Manual.mat'\n')
        return
    end
else
    if isfile('./ic_auto_selection/DecisionVariables_Auto.mat')
        load('./ic_auto_selection/DecisionVariables_Auto.mat')
    else
        fprintf('Cannot find ./ic_auto_selection/DecisionVariables_Auto.mat\n')
        return
    end
end

cleaned_dir = ['cleaned_', ic_select]; % based on auto vs manual cleaning

% if a filter is applied, smoothing should be s, also set if smothing is ns
if strcmp(filterchoice, 'bp')
    smoothchoice = 's_';
elseif strcmp(filterchoice, 'hp')
    smoothchoice = 's_';
else
    if strcmp(filterchoice, 's')
        smoothchoice = 's_';
    else
        smoothchoice = '';
    end
end

% set filter choice properly for file grabbing later. Default is bandpass.
if strcmp(filterchoice, 'bp')
    filterchoice = 'bp_';
    mean_signal_tag = 's_bp_';
    if ~exist('highfreq_cutoff', 'var')
        highfreq_cutoff = 0.15;
    end
    if ~exist('lowfreq_cutoff', 'var')
        lowfreq_cutoff = 0.008;
    end
elseif strcmp(filterchoice, 'hp')
    filterchoice = 'hp_';
    mean_signal_tag = 's_hp_';
    highfreq_cutoff = 1/(2*tr);
    if ~exist('lowfreq_cutoff', 'var')
        lowfreq_cutoff = 0.008;
    end
else
    filterchoice = '';
    mean_signal_tag = 's_hp_'; % easier to compare mean signals later if high pass filtered
    highfreq_cutoff = 1/(2*tr);
    lowfreq_cutoff = 0;
end

choice_tag = [smoothchoice, filterchoice];
% find/grab the correct comparison files - e.g., prefix_8p_funcfile_bold.nii.gz or
% prefix_9p_funcfile_bold.nii.gz
compare_filename_adj = [prefix, choice_tag, compare_file_tag, suffix]; % adj for adjusted by smoothing and filtering choices
compare_file = ['./', cleaned_dir, '/', compare_filename_adj];

if isfile(compare_file)
    meanplot_compare_file = ['./', cleaned_dir, '/', prefix, mean_signal_tag, compare_file_tag, suffix];
else
    compare_file = ['./cleaned_auto/', compare_filename_adj];
    if isfile(compare_file)
        meanplot_compare_file = ['./cleaned_auto/', prefix, mean_signal_tag, compare_file_tag, suffix];
    else
        compare_file = ['./', cleaned_dir, '/', compare_filename_adj];
        fprintf('Cannot find the compare file: \n')
        fprintf([compare_file, '\nExiting...\n'])
        return
    end
end



% grab original data for comparisons
orig_file_tag = 'orig';
orig_filename_adj = [prefix, smoothchoice, filterchoice, orig_file_tag, suffix];
orig_file = ['./', cleaned_dir, '/', orig_filename_adj];
meanplot_orig_file = ['./', cleaned_dir, '/', prefix, mean_signal_tag, orig_file_tag, suffix];

% Grab cadica data
denoised_tag = ['CADICA_', ic_select, '_', aggressionchoice];
denoised_filename = [prefix, denoised_tag, suffix];
denoised_filename_adj = [prefix, choice_tag, denoised_tag, suffix];
denoised_file = ['./', cleaned_dir, '/', denoised_filename_adj];
meanplot_denoised_file = ['./', cleaned_dir, '/', prefix, mean_signal_tag, denoised_tag, suffix];

% calculate final estimation of DOF - probably want > 15 DOF, and/or at
% least 10% of percent_variance_kept from CADICA
% something to keep in mind is that you should actually calculate the power
% percentage being removed to be much more accurate - need to find the
% power fraction of bp from the data that has been cleaned by CADICA
% without frequency filtering thus far. Would use denoised_filename without
% choice tag. Do later to improve this system.
Results.percent_freq_kept = ((highfreq_cutoff-lowfreq_cutoff) / (1/(2*tr)));
DOF_estimate_final = Data.numvolumes .* Results.percent_variance_kept .* Results.percent_freq_kept;
Results.num_ICs_kept = length(Results.signal_ICs);
Results.num_ICs_total = length(Results.ICs);

% confounds:
confound_place = './confounds_timeseries.csv';
allconfounds = readtable(confound_place);
confounds_dvars = table2array(allconfounds(:,{'dvars'}));
confounds_fd = table2array(allconfounds(:,{'framewise_displacement'}));
confounds_CSF = table2array(allconfounds(:,{'csf'}));
confounds_WM = table2array(allconfounds(:,{'white_matter'}));
confounds_GS = table2array(allconfounds(:,{'global_signal'}));

% Get GM prob mask and create GM mask and not GM mask for denoised
% and original files
fprintf(['Estimate of Final DOF for ', denoised_filename, ' is %.2f\n'], DOF_estimate_final)
GM_prob_file = './region_masks/GM_prop.nii.gz';
WMCSF_prob_file = './region_masks/WMandCSF_prop.nii.gz';
CSF_prob_file = './region_masks/CSF_prop.nii.gz';
Edge_prob_file = './region_masks/Edge_prop.nii.gz';
Outbrain_prob_file = './region_masks/OutbrainOnly_prop.nii.gz';
funcmask = './funcmask.nii.gz';
if isfile('./region_masks/Outbrain_prop.nii.gz') == 1
    NotGM_prob_file = './region_masks/Outbrain_prop.nii.gz'; % this is a better estimate if we calculated it
else
    NotGM_prob_file = './region_masks/NotGMorWM_prop.nii.gz'; % Makes sure not including anything that could be signal, removes WM too
end

% Now, get data for NotGM, Edge, Outbrain, WMCSF, and CSF for denoised data
fprintf('Calculating Denoised Data\n')
[denoised_GM, denoised_NotGM, denoised_GS, denoised_GM_mean, denoised_GM_diff_var, denoised_NotGM_diff_var] = getData(denoised_file, funcmask, GM_prob_file, NotGM_prob_file);
[~, denoised_Edge, ~, ~, ~, denoised_Edge_diff_var] = getData(denoised_file, funcmask, GM_prob_file, Edge_prob_file);
[~, denoised_Outbrain, ~, ~, ~, denoised_Outbrain_diff_var] = getData(denoised_file, funcmask, GM_prob_file, Outbrain_prob_file);
[~, denoised_WMCSF, ~, ~, ~, denoised_WMCSF_diff_var] = getData(denoised_file, funcmask, GM_prob_file, WMCSF_prob_file);
[~, denoised_CSF, ~, ~, ~, denoised_CSF_diff_var] = getData(denoised_file, funcmask, GM_prob_file, CSF_prob_file);

fprintf(['Running comparison of ', denoised_filename_adj, ' to ', compare_filename_adj, '\n'])    

% Now, get data for NotGM, Edge, Outbrain, WMCSF, and CSF for orig data
fprintf('Calculating Original Data\n')
[orig_GM, orig_NotGM, ~, ~, orig_GM_diff_var, orig_NotGM_diff_var] = getData(orig_file, funcmask, GM_prob_file, NotGM_prob_file);
[~, orig_Edge, ~, ~, ~, orig_Edge_diff_var] = getData(orig_file, funcmask, GM_prob_file, Edge_prob_file);
[~, orig_Outbrain, ~, ~, ~, orig_Outbrain_diff_var] = getData(orig_file, funcmask, GM_prob_file, Outbrain_prob_file);
[~, orig_WMCSF, ~, ~, ~, orig_WMCSF_diff_var] = getData(orig_file, funcmask, GM_prob_file, WMCSF_prob_file);
[~, orig_CSF, ~, ~, ~, orig_CSF_diff_var] = getData(orig_file, funcmask, GM_prob_file, CSF_prob_file);

% Now, get data for NotGM, Edge, Outbrain, WMCSF, and CSF for compared data
fprintf('Calculating Compare Data\n')
[compare_GM, compare_NotGM, ~, ~, compare_GM_diff_var, compare_NotGM_diff_var] = getData(compare_file, funcmask, GM_prob_file, NotGM_prob_file);
[~, compare_Edge, ~, ~, compare_Edge_diff_var] = getData(compare_file, funcmask, GM_prob_file, Edge_prob_file);
[~, compare_Outbrain, ~, ~, ~, compare_Edge_diff_var] = getData(compare_file, funcmask, GM_prob_file, Outbrain_prob_file);
[~, compare_WMCSF, ~, ~, ~, compare_Edge_diff_var] = getData(compare_file, funcmask, GM_prob_file, WMCSF_prob_file);
[~, compare_CSF, ~, ~, ~, compare_Edge_diff_var] = getData(compare_file, funcmask, GM_prob_file, CSF_prob_file);

% also get average GM signal of original and compare files (filtered for
% easier comparison)
fprintf('Calculating Mean GM Signals\n')
[~, ~, orig_GS_sbp, orig_GM_mean_sbp, ~, ~] = getData(meanplot_orig_file, funcmask, GM_prob_file, NotGM_prob_file);
[~, ~, compare_GS_sbp, compare_GM_mean_sbp, ~, ~] = getData(meanplot_compare_file, funcmask, GM_prob_file, NotGM_prob_file);
[~, ~, denoised_GS_sbp, denoised_GM_mean_sbp, ~, ~] = getData(meanplot_denoised_file, funcmask, GM_prob_file, NotGM_prob_file);

% Now compute relevant correlations, and some others that we don't use, but
% one could make use of if they desired
fprintf('Calculating Relevant Correlations\n')

% create permutations for regions, subtract a little bit because sometimes
% something funky happens and the sizes do not fully align with Edge
% Voxels. Very close (probably related to thresholding values)
GM_randperm = randperm(size(orig_GM, 1));
NotGM_randperm = randperm(size(orig_NotGM, 1));
Outbrain_randperm = randperm(size(orig_Outbrain, 1));
Edge_randperm = randperm(size(orig_Edge, 1));
WMCSF_randperm = randperm(size(orig_WMCSF, 1));
CSF_randperm = randperm(size(orig_CSF, 1));

% NotGM: Want to be shifted tighter and closer to 0
denoised_GM_NotGM_corr = createHistData(denoised_GM, denoised_NotGM, GM_randperm, NotGM_randperm, 500);
compare_GM_NotGM_corr = createHistData(compare_GM, compare_NotGM, GM_randperm, NotGM_randperm, 500);
orig_GM_NotGM_corr = createHistData(orig_GM, orig_NotGM, GM_randperm, NotGM_randperm, 500);

% GM: But, we want to maintain GM correlations in comparison (after other is removed)
% this one we compare to the original data, to see how much GM is
% maintained (want bad stuff removed, but not good stuff too!)
denoised_GM_GM_corr = createHistData(denoised_GM, orig_GM, GM_randperm, GM_randperm, 500);
compare_GM_GM_corr = createHistData(compare_GM, orig_GM, GM_randperm, GM_randperm, 500);
%orig_GM_GM_corr = createHistData(orig_GM, orig_GM, GM_randperm, GM_randperm, 500);

% Edge
denoised_GM_Edge_corr = createHistData(denoised_GM, denoised_Edge, GM_randperm, Edge_randperm, 500);
compare_GM_Edge_corr = createHistData(compare_GM, compare_Edge, GM_randperm, Edge_randperm, 500);
orig_GM_Edge_corr = createHistData(orig_GM, orig_Edge, GM_randperm, Edge_randperm, 500);

% Outbrain
denoised_GM_Outbrain_corr = createHistData(denoised_GM, denoised_Outbrain, GM_randperm, Outbrain_randperm, 500);
compare_GM_Outbrain_corr = createHistData(compare_GM, compare_Outbrain, GM_randperm, Outbrain_randperm, 500);
orig_GM_Outbrain_corr = createHistData(orig_GM, orig_Outbrain, GM_randperm, Outbrain_randperm, 500);

% WMCSF
denoised_GM_WMCSF_corr = createHistData(denoised_GM, denoised_WMCSF, GM_randperm, WMCSF_randperm, 500);
compare_GM_WMCSF_corr = createHistData(compare_GM, compare_WMCSF, GM_randperm, WMCSF_randperm, 500);
orig_GM_WMCSF_corr = createHistData(orig_GM, orig_WMCSF, GM_randperm, WMCSF_randperm, 500);

% CSF
denoised_GM_CSF_corr = createHistData(denoised_GM, denoised_CSF, GM_randperm, CSF_randperm, 500);
compare_GM_CSF_corr = createHistData(compare_GM, compare_CSF, GM_randperm, CSF_randperm, 500);
orig_GM_CSF_corr = createHistData(orig_GM, orig_CSF, GM_randperm, CSF_randperm, 500);

% dvars
denoised_GM_dvars_corr = createHistConf1DData(denoised_GM, confounds_dvars, GM_randperm, 10000);
compare_GM_dvars_corr = createHistConf1DData(compare_GM, confounds_dvars, GM_randperm, 10000);
orig_GM_dvars_corr = createHistConf1DData(orig_GM, confounds_dvars, GM_randperm, 10000);

% fd
denoised_GM_fd_corr = createHistConf1DData(denoised_GM, confounds_fd, GM_randperm, 10000);
compare_GM_fd_corr = createHistConf1DData(compare_GM, confounds_fd, GM_randperm, 10000);
orig_GM_fd_corr = createHistConf1DData(orig_GM, confounds_fd, GM_randperm, 10000);


% idea take voxels, run correlation of voxels to themselves, but before and
% after denoising. NotGM should be less correlated than GM
% denoised_GM_autovoxel_corr = autovoxelcorr(denoised_GM, orig_GM, GM_randperm, 1000);
% denoised_NotGM_autovoxel_corr = autovoxelcorr(denoised_NotGM, orig_NotGM, NotGM_randperm, 1000);
% compare_GM_autovoxel_corr = autovoxelcorr(compare_GM, orig_GM, GM_randperm, 1000);
% compare_NotGM_autovoxel_corr = autovoxelcorr(compare_NotGM, orig_NotGM, NotGM_randperm, 1000);



%s_comparison_string = [s_compare_filename, ' vs denoised'];
%s_filterchoice_comparison_string = [s_filterchoice_compare_filename, ' vs denoised'];

% Create QC folder (if it doesn't already exist)
if ~isfolder('qc')
    mkdir qc
end
cd qc

% Save relevant QC Data so you do not have to redo this later
% parse compare filename first 
splitting_compare = split(compare_file_tag, '_');

% if the compare string already has subject information in it (like auto
% CADICA, we do not need that information

parsed_compare = join(splitting_compare(length(split(prefix, '_')):end-(length(split(suffix, '_'))-1)), '_'); % gets rid of prefix and suffix tags for each file
parsed_compare = parsed_compare{:};

split_prefix = split(prefix, '_');
prefix_cut = strjoin(split_prefix(1:end-1), '_');
split_choice_tag = split(choice_tag, '_');
choice_tag_cut = strjoin(split_choice_tag(1:end-1), '_'); % gets rid of final underscore


qc_naming = [prefix, choice_tag, denoised_tag, '_vs_', compare_file_tag]; 
%qc_naming = [smoothchoice, filterchoice, ic_select, '_', aggressionchoice, '_vs_', parsed_compare];
qc_vals = [qc_naming, '_qc_vals.mat'];
qc_plots = [qc_naming, '_qc_plots.png'];

% save relevant plotting variables
fprintf('Saving Relevant QC Data\n')
title_string = [prefix_cut, ': CADICA (', smoothchoice, filterchoice, ic_select, '_', aggressionchoice, '), ', compare_file_tag, ', & ' smoothchoice, filterchoice,'orig'];
%save(qc_vals)
save(qc_vals, 'orig_GM_NotGM_corr', 'compare_GM_NotGM_corr', 'denoised_GM_NotGM_corr', ...
     'compare_GM_GM_corr', 'denoised_GM_GM_corr', ...
     'orig_GM_Outbrain_corr', 'compare_GM_Outbrain_corr', 'denoised_GM_Outbrain_corr', ...
     'orig_GM_Edge_corr', 'compare_GM_Edge_corr', 'denoised_GM_Edge_corr', ...
     'orig_GM_WMCSF_corr', 'compare_GM_WMCSF_corr', 'denoised_GM_WMCSF_corr', ...
     'orig_GM_CSF_corr', 'compare_GM_CSF_corr', 'denoised_GM_CSF_corr', ...
     'orig_GM_dvars_corr', 'compare_GM_dvars_corr', 'denoised_GM_dvars_corr', ...
     'orig_GM_fd_corr', 'compare_GM_fd_corr', 'denoised_GM_fd_corr', ...
     'orig_GM_mean_sbp', 'compare_GM_mean_sbp', 'denoised_GM_mean', ...
     'title_string', 'orig_file_tag', 'compare_file_tag', 'denoised_filename',...
     'filterchoice', 'aggressionchoice', 'ic_select', 'smoothchoice', 'mean_signal_tag', ...
     'DOF_estimate_final', 'Results', 'tr')

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
histogram(orig_GM_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_fd corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against overall movement.
nexttile
hold on
title('FD-GM Corr ', 'Interpreter', 'none')
histogram(orig_GM_fd_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_fd_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_fd_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_dvars corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against global signal changes (including large
% movement)
nexttile
hold on
title('DVARS-GM Corr ', 'Interpreter', 'none')
histogram(orig_GM_dvars_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_dvars_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_dvars_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_Outbrain corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('Outbrain-GM Corr ', 'Interpreter', 'none')
histogram(orig_GM_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_WMCSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('WMCSF-GM Corr ', 'Interpreter', 'none')
histogram(orig_GM_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_CSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('CSF-GM Corr ', 'Interpreter', 'none')
histogram(orig_GM_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_NotGM corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that. Want to see the denoised and orig be significantly different
nexttile
hold on
title('NotGM-GM Corr ', 'Interpreter', 'none')
histogram(orig_GM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(compare_GM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
histogram(denoised_GM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
hold off

% GM_GM autocorr to orig: Want to see signal correlations maintained, in comparison to
% GM_NotGM corr, but still reduced (removed erroneous/global correlations).
% get color order right to match other plots
colorord = get(gca, 'colororder');

nexttile
hold on
title('GM-GM Auto Corr to Orig', 'Interpreter', 'none')
%histogram(orig_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs')
histogram(compare_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
histogram(denoised_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
legend('compare', 'cadica', 'Interpreter', 'none')
hold off

% Plot GS changes before and after Denoising - it is easier to
% compare a smoothed and bandpassed raw data to a smooth and
% bandpassed cadica data (easier to tell if/how major noise
% spikes may be removed by cadica). Suggests noise spiking is
% likely greatly reduced, and helps show if GS signal is
% reasonable (did something go horribly wrong?)
nexttile
hold on
title('Mean GM Signal (Filtered)', 'Interpreter', 'none')
plot(orig_GM_mean_sbp, 'LineWidth', 1.5)
plot(compare_GM_mean_sbp, 'LineWidth', 1.5)
plot(denoised_GM_mean, 'LineWidth', 1.5)
legend('orig', 'compare', 'cadica', 'Interpreter', 'none')
xlim([0, length(orig_GM_mean_sbp)])
hold off


fprintf('Saving Figure\n')
exportgraphics(t, qc_plots, 'Resolution', 600)

fprintf('\n')
end

%%
function compare_autovoxel_corr = autovoxelcorr(compare_region, orig_region, region_randperm, perms)
% Run correlations of voxels to themselves from before to after denoising. 
% Purpose: Good denoising should maintain more GM autocorrelation than
% notGM, in theory
    compare_autovoxel_corr = corr(compare_region(region_randperm(1:perms), :)', orig_region(region_randperm(1:perms), :)');
    compare_autovoxel_corr = compare_autovoxel_corr(tril(compare_autovoxel_corr, -1) ~= 0);
end

function [GM, region, Global_Signal, GM_Signal_mean, GM_var, region_var] = getData(funcfile, funcmask, GM_prob_file, region_prob_file)
% read in the niftis
region_prob_data = niftiread(region_prob_file);
gm_prob_data = niftiread(GM_prob_file);
func_data = niftiread(funcfile);
funcmask_data = niftiread(funcmask);

gm_mask = logical(gm_prob_data > 0.75 & funcmask_data == 1);
region_mask = logical(region_prob_data > 0.75 & funcmask_data == 1);
func_mask = logical(funcmask_data == 1);

% Select the timeseries for GM and region and global signal
global_signal = zeros(sum(func_mask(:)), size(func_data, 4));
GM = zeros(sum(gm_mask(:)), size(func_data, 4));
region = zeros(sum(region_mask(:)), size(func_data, 4));
for j = 1:size(func_data, 4)
    curr_func = func_data(:,:,:,j);
    GM(:,j) = curr_func(gm_mask);
    region(:,j) = curr_func(region_mask);
    global_signal(:,j) = curr_func(func_mask);
end

GM = GM(GM(:,1) ~= 0, :);
region = region(region(:, 1) ~= 0, :);
Global_Signal = global_signal(global_signal(:, 1) ~= 0, :);

GM_var = var(GM, 0, 2); % variance of difference in time for GM)
region_var = var(region, 0, 2);

Global_Signal_mean = mean(global_signal(global_signal(:, 1) ~= 0, :));
Global_Signal_mean = Global_Signal_mean - mean(Global_Signal_mean);
GM_Signal_mean = mean(GM(GM(:,1) ~= 0, :));
GM_Signal_mean = GM_Signal_mean - mean(GM_Signal_mean);
end

function compare_orig_corr = createHistData(compare_GM, orig_region, GM_randperm, region_randperm, perms)
% Should center close to 0 (not globally correlated) if it's
% denoised (e.g., denoised lowered noise)
compare_orig_corr = corr(compare_GM(GM_randperm(1:perms), :)', orig_region(region_randperm(1:perms), :)');
compare_orig_corr = compare_orig_corr(tril(compare_orig_corr, -1) ~= 0);

% Should center above 0 (e.g., denoise did not remove true signal)
% GM_corr = corr(compare_GM(GM_randperm(1:perms), :)', compare_GM(GM_randperm(1:perms), :)');
% GM_GM_corr = GM_corr(tril(GM_corr, -1) ~= 0);
end

function [GM_conf1D_corr] = createHistConf1DData(GM, conf, GM_randperm, perms)
% This is for confounds that are calculated in the 1st derivative
% and absolute value (e.g., dvars, fd)
% corr conf vs GM diff abs: Should be reduced if denoised
% consider detrend() for GM, since that is what we do for script 2 to focus
% more on spikes, less on drift.
GM_conf1D_corr = corr(conf(2:end), abs(diff(GM(GM_randperm(1:perms), :)')))';
end

function [GM_conf_corr] = createHistConfData(GM, conf, perms)
% This is for confounds that are calculated without derivative
% (e.g., CSF, WM, GS)
GM_randperm = randperm(size(GM, 1));
% corr conf vs GM: Should be reduced if denoised
GM_conf_corr = corr(conf, GM(GM_randperm(1:perms), :)')';
end