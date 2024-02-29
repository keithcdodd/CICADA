function cicada_group_qc(output_dir, cicada_dir, cicada_type, subj_ids, ses_ids, task_name, bad_data_list, manual_adjustment_list)
% function to run group qc on cicada output
% output_dir: where you want the group_qc: commonly something like: /cicada_group_qc/task_name
% cicada_dir: general cicada directory to read from
% cicada_type: 'auto' or 'manual'
% subj_ids: e.g., {'102', '102', '103', '103'}
% ses_ids: e.g., {'01', '02', '01', '02'}
% task_name (include run in there if it exists) e.g., 'visual_run-01'
% Note: Only one task_name 
% bad_data_list: (clearly something wrong with scan, unusable data, found
% earlier, like with mriqc). '1' designates bad. '0' is fine: e.g., {'0', '0', '1', '0'}
% manual_adjustment_list: Data you have tried to save with manual
% adjustment of IC selection (and Manual CICADA has been run): e.g., {'1',
% '0'. '1', '0'}
% Note: subj_ids, ses_ids, bad_data_list, and manual_adjustment_list must be all the same length of cell
% arrays! Easier to read this from a .csv you make and feed that into this
% script


cicada_group_qc_dir = fileparts(mfilename('fullpath')); % this gives current script path
cd([cicada_group_qc_dir, '/..'])
cicada_path = pwd;
gm_mni_prob_file = [cicada_path, '/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_gm_tal_nlin_asym_09c_2mm.nii.gz'];
background_file = [cicada_path, '/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2mm.nii.gz'];
network_file = [cicada_dir, '/templates/network_template_', task_id, '.nii.gz'];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% initialize struct:
Group_QC = struct;

% convert gm probability into a mask at .67
gm_mni_prob = niftiread(gm_mni_prob_file);
gm_mni_thresh = gm_mni_prob_file > 0.67;
gm_mni_thresh_info = niftiinfo(gm_mni_prob_file);
gm_mni_thresh_info.Datatype = 'uint8'; % in case I want to write it out as a file later
compare_tags = {'8p'};

% find the correct cicada cleaned folder tag, default is manual
if strcmp(cicada_type, 'auto') ~= 1
    cicada_type = 'manual';
    compare_tags = {'8p', 'auto'}; % if doing all manual, it is good to compare to all auto to see impact
end
cicada_type_default = cicada_type; % for converting later

% Create output_dir, if it does not already exist:
if not(isfolder(output_dir))
    mkdir(output_dir)
end


data_fol = [output_dir, '/data'];
if not(isfolder(data_fol))
    mkdir(data_fol)
else
   rmdir(data_fol, 's') % remove data folder and contents if it already exists
   mkdir(data_fol) % and now make new folder
end

% Make separate photo folders, so it is easier to just scroll through them
% all in the future
for h = 1:length(compare_tags)
    % create folders for comparison photos to zip through too
    photo_fol = [output_dir, '/qc_photos_cicada_', cicada_type, '_', compare_tags{h}];
    if not(isfolder(photo_fol))
        mkdir(photo_fol)
    else
        rmdir(photo_fol, 's')
        mkdir(photo_fol)
    end
end

% consider initializing a 4D signaltonoise file with zeros  to save on
% time
m = 1; % counting purposes, nothing major or wild

for j = 1:length(subj_ids)

    fprintf('\n')
    fprintf(['Running for subject ', subj_ids{j}, '\n'])
    % check for folder
    curr_subj_fol = [cicada_dir, '/sub-', subj_ids{j}];
    if ~isfolder(curr_subj_fol)
        fprintf(['No subject folder found for ', subj_ids{j}, '. Moving on to next one...\n'])
        continue
    end

    for k = 1:length(ses_ids)
        
        fprintf(['Running for session ', ses_ids{k}, '\n'])
        % check for folder
        curr_sess_fol = [curr_subj_fol, '/ses-', ses_ids{k}];
        if ~isfolder(curr_sess_fol)
            fprintf(['No sess folder found for ', ses_ids{k}, '. Moving on to next one...\n'])
            continue
        end
        
        % Can finally now get to the task
        fprintf(['Running for task ', task_id, '\n'])
        curr_task_fol = [curr_sess_fol, '/', task_id];

        if ~isfolder(curr_task_fol)
            fprintf(['No task folder found for ', task_id, '. Moving on to next one...\n'])
            continue
        end

        % OK, now grab the things you need/want:
        % The actually cleaned file, in cleaned_manual folder:
        prefix = ['sub-', subj_ids{j}, '_ses-', ses_ids{k}, '_task-', task_id];
        suffix = 'bold.nii.gz';

        % if it is labeled as bad data from the get-go (e.g., mriqc shows
        % that the data was corrupted from the start), don't include in
        % analysis
        % This if statement makes sure it is NOT the terrible data 
        if sum(contains(bad_data_prefixes, prefix)) < 1
            % cd into task folder and go from there:
            cd(curr_task_fol)

            % check if this is a subject where you want to use the manually
            % adjusted, or not include it
            % To Do: figure out what you want to do here then
            adjusted = 0; % 1 if we adjusted it
            cicada_type = cicada_type_default; % update to default version (often this is auto)
            denoise_tag = ['cicada_', cicada_type, '_', reg_type];
            % if the current subject falls under the list of adjust or bad
            % prefixes, grab the manual version and label it appropriately
            if sum(contains(adjust_prefixes, prefix)) > 0
                cicada_type = 'manual';
                denoise_tag = ['cicada_', cicada_type, '_', reg_type];
                adjusted = 1;
            end
    
             % check for cleaned folder
            cleaned_folname = ['cleaned_' cicada_type];
            cleaned_fol = [curr_task_fol, '/', cleaned_folname];
            if not(isfolder(cleaned_fol))
                fprintf(['No cleaned folder at ', cleaned_fol, ' found. Moving on to next one...\n'])
                continue
            end
    
            curr_filename = [prefix, '_', smooth_filter_tag, denoise_tag, '_', suffix];
            curr_denoised_file = [cleaned_fol, '/', curr_filename];
           
    
            % the funcmask
            curr_mask_file = [curr_task_fol, '/funcmask.nii.gz'];
    
            if ~exist(curr_denoised_file, 'file')
                fprintf(['No file matching ' curr_filename, ' in cleaned folder ' cleaned_fol, ' found. Moving on to next one...\n'])
            end
    
            % OK, we have the denoised file! Now copy it to group folder!
            copyfile(curr_denoised_file, data_fol)
            curr_copied_denoised_file = [data_fol, '/' curr_filename];
    
            % Grab current functional mask and load it into a larger 4D
            % file
            curr_funcmask = niftiread(curr_mask_file);
            funcmask(:,:,:,m) = curr_funcmask;
    
            % Grab Signal IC Overlap, found in ic_manual_selection, or
            % ic_auto_selection depending on your selection. Help you
            % determine if regions of interest are well represented per
            % subject
            curr_signalandnoise_overlap_filename = 'SignalandNoiseICOverlap.nii.gz'; % 1 for signal, -1 for noise without signal
    
            %curr_signalandnoise_overlap_filename = 'SignalICOverlap.nii.gz';
            
            selection_fol = [curr_task_fol, '/ic_', cicada_type, '_selection'];
            curr_signalandnoise_file = [selection_fol, '/', curr_signalandnoise_overlap_filename];
    
            curr_signalandnoise_mask = niftiread(curr_signalandnoise_file);
    
            % Calculate overlap with GM mask
            if m == 1
                flirt_command = ['flirt -ref ', curr_signalandnoise_file, ' -in ', gm_mni_prob_file, ' -out ', output_dir, '/gm_mni_prob_resampled.nii.gz -usesqform -applyxfm'];
                [status, cmdout_flirt] = system(flirt_command, '-echo');
                gm_mni_prob_resampled = niftiread([output_dir, '/gm_mni_prob_resampled.nii.gz']);
                gm_mni_thresh = gm_mni_prob_resampled > 0.67;
            end
            
            signal_mask = curr_signalandnoise_mask == 1;
            gm_signal_overlap = gm_mni_thresh .* signal_mask;
            
            gm_signal_proportion = sum(gm_signal_overlap(:)) / sum(gm_mni_thresh(:)); % What proportion of gm is covered by signal
            signal_gm_proportion = sum(gm_signal_overlap(:)) / sum(signal_mask(:)); % what proportion of signal is within GM
    
            % easier yet, calculate a dice coefficient
            gm_and_signal = gm_signal_overlap;
            gm_or_signal = (gm_mni_thresh + signal_mask) > 0;
            gm_signal_dice = sum(gm_and_signal(:)) / sum(gm_or_signal(:));
    
            % load it into larger 4D file
            signalandnoise_mask(:,:,:,m) = curr_signalandnoise_mask;
    
            % load the qc values, and then save the ones you care about
            % most
            qc_fol = [curr_task_fol, '/qc'];
            qc_vals = [qc_fol, '/', prefix, '_', denoise_tag, '_vs_', compare_tags{1}, '_qc_vals.mat']; % load values of interest will be the same independent of compare tag, so just use first one
            load(qc_vals)
            percent_ICs_kept = Results.num_ICs_kept / Results.num_ICs_total;
    
            % OK we definitely want to save the following in a table(things that are one value per task): num_ICs_total,
            % num_ICs_kept, DOF_estimate_final
            % if no filtering is applied, dof_estimate_final will be equal to percent_variance_kept
            
            % dof_estimate_final should be equal to percent_variance_kept *
            % percent_freq_kept * numvolumes. Percent frequency kept is
            % 100% if we did not do any frequency filtering
    
            % in case we did not do "hisorically standard" resting state bp
            % (0.01-0.1Hz), calculate what estimated DOF would become if we
            % did do it.
            percent_freq_kept_if_standardbp = 0.09 / (1 / (2* tr));
            DOF_estimate_if_standard_bp = Results.percent_variance_kept * percent_freq_kept_if_standardbp * numvolumes;
    
            % Get confounds too:
            confound_place = [curr_task_fol, '/confounds_timeseries.csv'];
            allconfounds = readtable(confound_place);
            FD = table2array(allconfounds(:,{'framewise_displacement'}));
            DVARS = table2array(allconfounds(:,{'dvars'}));
            RMS = table2array(allconfounds(:,{'rmsd'})); 
    
            medFD = median(FD(2:end)); % median FD to pick cut offs. Liberal is 0.55mm cut off, conservative is >0.25mm cut off
            meanFD = mean(FD(2:end)); % some QC_FC likes comparing to mean FD
            perc_FD_above_thresh = (sum(FD(2:end) > 0.2) / length(FD(2:end))) * 100; %percent FD greater 0.2mm, 20% cut off also used for conservative
            Any_FD_above_thresh = int8(sum(FD(2:end) > 5) > 0); % also added onto conservative cut off, any FD > 5mm
            medDVARS = median(FD(2:end));
            meanRMS = mean(RMS(2:end)); % Add this into qc_group_array
    
            % And lets include mean and std for all relevant QC plots from script 4:
            GM_NotGM_mean = median(denoised_GM_NotGM_corr);
            GM_NotGM_sd = std(denoised_GM_NotGM_corr);
            GM_GM_mean = median(denoised_GM_GM_corr);
            GM_GM_sd = std(denoised_GM_GM_corr);
            GM_dvars_mean = median(denoised_GM_dvars_corr);
            GM_dvars_sd = std(denoised_GM_dvars_corr);
            GM_fd_mean = median(denoised_GM_fd_corr);
            GM_fd_sd = std(denoised_GM_fd_corr);
            GM_CSF_mean = median(denoised_GM_CSF_corr);
            GM_CSF_sd = std(denoised_GM_CSF_corr);
            GM_WMCSF_mean = median(denoised_GM_WMCSF_corr);
            GM_WMCSF_sd = std(denoised_GM_WMCSF_corr);
            GM_Outbrain_mean = median(denoised_GM_Outbrain_corr);
            GM_Outbrain_sd = std(denoised_GM_Outbrain_corr);
            GM_Edge_mean = median(denoised_GM_Edge_corr);
            GM_Edge_sd = std(denoised_GM_Edge_corr);
    
            qc_group_array(:,m) = [{m, curr_filename, subj_ids{j}, ses_ids{k}, task_id, adjusted, meanRMS, ...
                medFD, meanFD, perc_FD_above_thresh, Any_FD_above_thresh, medDVARS, gm_signal_dice, gm_signal_proportion, Results.num_ICs_kept, ...
                Results.num_ICs_total, percent_ICs_kept, Results.percent_variance_kept, ...
                Results.percent_freq_kept, numvolumes, DOF_estimate_final, DOF_estimate_if_standard_bp, ...
                GM_NotGM_mean, GM_NotGM_sd, GM_GM_mean, GM_GM_sd, GM_dvars_mean, GM_dvars_sd, ...
                GM_fd_mean, GM_fd_sd, GM_CSF_mean, GM_CSF_sd, GM_WMCSF_mean, GM_WMCSF_sd, ...
                GM_Outbrain_mean, GM_Outbrain_sd, GM_Edge_mean, GM_Edge_sd}, num2cell(Results.compare_cleaning.After')];
    
            data_files{m} = curr_copied_denoised_file;
    
            % finally copy over the photos to their respective folders
            for h = 1:length(compare_tags)
                photo_fol = [output_dir, '/qc_photos_', compare_tags{h}];
                photo_name = [prefix, '_', denoise_tag, '_vs_', compare_tags{h}, '_qc_plots.png'];
                photo_file = [curr_task_fol, '/qc/', photo_name];
                copyfile(photo_file, photo_fol)
            end
            m = m+1; % increment m
        end

    end
end

% record the bad data that was not even looked at:
Group_QC.bad_data_prefixes = bad_data_prefixes;

% go to general group folder for this specific set
cd(output_dir)

% cell to table
qc_table_labels = ['image_number', 'image_names', 'subject', 'session', 'task', 'manually_adjusted', 'meanRMS', ...
    'median_FD', 'mean_FD', 'Percent_FD_gt_point2mm','AnyFD_gt_5mm', 'median_DVARS', 'gm_signal_dice', 'gm_coverage', 'number_kept_ics', 'number_total_ics', ...
    'fraction_kept_ics', 'fraction_signal_variance_kept', 'fraction_frequency_kept', 'numvolumes', ...
    'dof_estimate_final', 'dof_estimate_final_if_standard_resting_bp', 'GM_NotGM_median', 'GM_NotGM_sd', ...
    'GM_GM_median', 'GM_GM_sd', 'GM_dvars_median', 'GM_dvars_sd', 'GM_fd_median', 'GM_fd_sd', ...
    'GM_CSF_median', 'GM_CSF_sd', 'GM_WMCSF_median', 'GM_WMCSF_sd', 'GM_Outbrain_median', 'GM_Outbrain_sd', ...
    'GM_Edge_median', 'GM_Edge_sd', Results.compare_cleaning.Properties.RowNames'];
final_qc_table = cell2table(qc_group_array', 'VariableNames', qc_table_labels);

% Calculate outliers based on the final qc table. Do it in three ways and
% add it to the qc table:

% Number of total ICs (if very low -- bad data all round)
Group_QC.low_number_total_ics = (isoutlier(final_qc_table.number_total_ics, "median")) & (final_qc_table.number_total_ics < mean(final_qc_table.number_kept_ics));

% Low GM coverage, regular and dice
Group_QC.low_gm_coverage = (isoutlier(final_qc_table.gm_coverage, "median")) & (final_qc_table.gm_coverage < mean(final_qc_table.gm_coverage));
Group_QC.low_gm_dice = (isoutlier(final_qc_table.gm_signal_dice, "median")) & (final_qc_table.gm_signal_dice < mean(final_qc_table.gm_signal_dice));

% Fraction Signal Variance Kept (low would suggest there are very few good 
% number of ICs (too swamped by noise) and it does not explain much of the data):
Group_QC.low_fraction_signal_variance_kept = (isoutlier(final_qc_table.fraction_signal_variance_kept, "median")) & (final_qc_table.fraction_signal_variance_kept < mean(final_qc_table.fraction_signal_variance_kept));

% Signal (The kept ICs should contain a fair amount of signal region
% overlap)
Group_QC.low_Signal = (isoutlier(final_qc_table.Signal, "median")) & (final_qc_table.Signal < mean(final_qc_table.Signal));

% Smoothing
Group_QC.low_Smoothing = (isoutlier(final_qc_table.Smoothing_Retention, "median")) & (final_qc_table.Smoothing_Retention < mean(final_qc_table.Smoothing_Retention));

% Power Overlap
% if task data, take into account task power overlap
Group_QC.low_general_power_overlap = (isoutlier(final_qc_table.general_power_overlap, "median")) & ...
    (final_qc_table.general_power_overlap < mean(final_qc_table.general_power_overlap));
if has_task_onsets == 1
    Group_QC.low_besttask_power_overlap = (isoutlier(final_qc_table.besttask_power_overlap, "median")) & ...
    (final_qc_table.besttask_power_overlap < mean(final_qc_table.besttask_power_overlap));
    Group_QC.low_power_overlap = Group_QC.low_general_power_overlap & Group_QC.low_besttask_power_overlap;
else
    %else it is resting state
    Group_QC.low_power_overlap = Group_QC.low_general_power_overlap;
end


% DVARS_Corr
Group_QC.high_DVARS_corr = (isoutlier(abs(final_qc_table.GM_dvars_median), "median")) & (abs(final_qc_table.GM_dvars_median) > median(abs(final_qc_table.GM_dvars_median)));

% FD_Corr
Group_QC.high_FD_corr = (isoutlier(abs(final_qc_table.GM_fd_median), "median")) & (abs(final_qc_table.GM_fd_median) > median(abs(final_qc_table.GM_fd_median)));


% combine to get overall outliers
Group_QC.cicada_outliers = logical(Group_QC.low_number_total_ics + ...
    Group_QC.low_fraction_signal_variance_kept + Group_QC.low_Signal + ...
    Group_QC.low_Smoothing + Group_QC.low_power_overlap + Group_QC.low_gm_coverage + Group_QC.low_gm_dice);

% Now, get liberal outliers (not restrictive)
Group_QC.liberal_outliers = final_qc_table.mean_FD > 0.55;

% Finally, conservative outliers: median FD
Group_QC.conservative_outliers = (final_qc_table.mean_FD > 0.25) | ...
    (final_qc_table.Percent_FD_gt_point2mm > 20) | ...
    (final_qc_table.AnyFD_gt_5mm);

final_qc_table.low_number_total_ics = Group_QC.low_number_total_ics;
final_qc_table.low_fraction_signal_variance_kept = Group_QC.low_fraction_signal_variance_kept;
final_qc_table.low_Signal = Group_QC.low_Signal;
final_qc_table.low_Smoothing = Group_QC.low_Smoothing;
final_qc_table.low_power_overlap = Group_QC.low_power_overlap;
final_qc_table.low_gm_coverage = Group_QC.low_gm_coverage;
final_qc_table.high_DVARS_corr = Group_QC.high_DVARS_corr;
final_qc_table.high_FD_corr = Group_QC.high_FD_corr;
final_qc_table.cicada_outliers = Group_QC.cicada_outliers;
final_qc_table.liberal_outliers = Group_QC.liberal_outliers;
final_qc_table.conservative_outliers = Group_QC.conservative_outliers;

writetable(final_qc_table, 'group_qc_table.csv') % can use this and sort by gm_signal and signal_gm proportions and also number_kept_ics to help find the bad or adjustable data
Group_QC.final_qc_table = final_qc_table;

% Now work on the images:

% remove outliers from funcmask 
image_keep = ~Group_QC.cicada_outliers;

% remove outliers from signal and noise mask
signalandnoise_mask_all = signalandnoise_mask;
signalandnoise_mask = signalandnoise_mask_all(:,:,:,image_keep);

% collapse funcmask into 3D with mean, remove outliers
mean_funcmask = mean(funcmask(:,:,:,image_keep), 4); 
funcmask_info = niftiinfo(curr_mask_file); % will be 3D instead of 4D, will need 
% thresh is for melodic mask - mask 99% coverage needed
funcmask_info.Datatype = 'uint8';
niftiwrite(cast((mean_funcmask>0.99), 'uint8'), 'mean_funcmask_thresh', funcmask_info, 'Compressed', true) % 0-1, it should be the percent of funcmasks that overlapped with voxel
% prob is for helping with a final group mask
funcmask_info.Datatype = 'single';
niftiwrite(cast(mean_funcmask, 'single'), 'mean_funcmask_prob', funcmask_info, 'Compressed', true) % 0-1, it should be the percent of funcmasks that overlapped with voxel

% Create 4d funcmask too so it is easy to search through
funcmask_info.Datatype = 'uint8';
funcmask_info.ImageSize = [funcmask_info.ImageSize, size(signalandnoise_mask,4)];
funcmask_info.PixelDimensions = [funcmask_info.PixelDimensions, 1];
niftiwrite(cast(funcmask(:,:,:,image_keep), 'uint8'), 'funcmask_Group', funcmask_info, 'Compressed', true)

% OK, now load and create 4D signaltonoise file and header, and write to
% the group directory
signal_prob_info = niftiinfo(curr_signalandnoise_file); % will be 3D instead of 4D, will need to edit to 4D
signal_prob_info.ImageSize = [signal_prob_info.ImageSize, size(signalandnoise_mask,4)];
signal_prob_info.PixelDimensions = [signal_prob_info.PixelDimensions, 1]; % idk if 1 is "right", but it makes some sense/should work! 4 dimension is just each image

niftiwrite(cast(signalandnoise_mask, 'single'), 'SignalandNoise_Group', signal_prob_info, 'Compressed', true) % 1 is high signal, -1 is high noise-labeled only
%niftiwrite(cast(signalandnoise_mask, 'single'), 'Signal_Group', signal_prob_info, 'Compressed', true) % 1 is high signal, -1 is high noise-labeled only



% Run Group MELODIC
% First, make text file of the image names
% will need a text file where separated by commas (no spaces) and then in
% the melodic call $(cat )
% Make sure to not include data files that were designated as bad enough to
% remove
all_data_files = data_files;
data_files = all_data_files(~Group_QC.cicada_outliers);
writecell(data_files, 'image_names.txt')

% Resample mni to be the background image for melodic
flirt_command = ['flirt -ref mean_funcmask_thresh.nii.gz -in ', background_file, ' -out mni_resampled.nii.gz -usesqform -applyxfm'];
[status, cmdout_flirt] = system(flirt_command, '-echo');

% If you want to specify how many ICs to make:
%Group_ICA_command = ['melodic --in=$(cat image_names.txt) --outdir=melodic --Ostats -m mean_funcmask_thresh.nii.gz --bgimage=mni_resampled.nii.gz -d ', num2str(num_ICs), ' --nobet --mmthresh=0.5 --report --tr=', num2str(tr)];
% Vs.: let Melodic decide on number of ICs (preferred)
Group_ICA_command = ['melodic --in=$(cat image_names.txt) --outdir=melodic --Ostats -m mean_funcmask_thresh.nii.gz --bgimage=mni_resampled.nii.gz --nobet --mmthresh=0.5 --report --tr=', num2str(tr)];

if redo_melodic == 0
    if not(isfolder('melodic'))
        fprintf(['Running: ', Group_ICA_command, '\n'])
        %[status, cmdout_ICA] = system(Group_ICA_command, '-echo');
        fprintf('Done with Group ICA (Melodic)\n\n')
    end
elseif redo_melodic == 1
    fprintf(['Running: ', Group_ICA_command, '\n'])
    [status, cmdout_ICA] = system(Group_ICA_command, '-echo');
    fprintf('Done with Group ICA (Melodic)\n\n')

    % create ICprobabilities combined in melodic folder (good for making masks
    % later at 100% probability, for example).
    % first get to qc folder in bash
    out = system(['cd ', output_dir], '-echo');
    [~, str] = system('ls ./melodic/stats/probmap_* | sort -V');
    probmapnames = regexprep(str, '\s+', ' '); % convert newlines to spaces
    command = ['fslmerge -t ./melodic/ICprobabilities.nii.gz ', probmapnames];
    [status, cmdout] = system(command, '-echo');
    % And load this now, 4th dimension is the ICs
    IC_probs = niftiread('./melodic/ICprobabilities.nii.gz');
    IC_probs_info = niftiinfo('./melodic/ICprobabilities.nii.gz');
    IC_probs_99 = logical(IC_probs > 0.99); % keeps only very best overlap
    
    % Update headers for 3D images
    IC_3D_info = IC_probs_info;
    IC_3D_info.ImageSize = IC_3D_info.ImageSize(1:end-1);
    IC_3D_info.PixelDimensions = IC_3D_info.PixelDimensions(1:end-1);
    IC_3D_info.Datatype = 'single';
    
    
    % Also go through melodic output and calculate overlap % of IC to each
    % network between the 7 networks, which can be useful
    networks_combined = niftiread(network_file);
    networks_info = niftiinfo(network_file);
    
    % Break it up into the 7 networks in the 4th dimension
    networks = zeros([size(networks_combined),7]);
    for i1 = 1:7
        networks(:,:,:,i1) = logical(networks_combined == i1);
    end
    
    % OK, now you can just multiply this by ICprobabilities file since each 4th
    % dimension is an IC
    % because matlab does not love a 4th dimension, the code is a little
    % complex
    % initialize:
    IC_in_network = zeros(size(IC_probs_99,4), size(networks,4));
    network_in_IC = IC_in_network;
    dice_IC_network = IC_in_network;
    for i1 = 1:size(networks,4)
        curr_network = networks(:,:,:,i1);
    
        curr_IC_network_overlap = zeros(size(IC_probs_99));
        for i2 = 1:size(IC_probs_99,4)
            % make matrix only showing overlap between probability maps and
            % current network
            curr_IC_network_overlap(:,:,:,i2) = logical(IC_probs_99(:,:,:,i2) .* curr_network); % Try to mask in 4th dimension
        end
        
        % sum across everything except ICs & squeeze to reduce dimensions to
        % calculate proportional overlap per IC
        IC_in_network(:,i1) = squeeze(sum(sum(sum(curr_IC_network_overlap,1),2),3)) ./ squeeze(sum(sum(sum(IC_probs_99,1),2),3));
        network_in_IC(:,i1) = squeeze(sum(sum(sum(curr_IC_network_overlap,1),2),3)) / sum(curr_network(:));
    
        dice_IC_network(:,i1) = squeeze(sum(sum(sum(curr_IC_network_overlap,1),2),3)) ./ ...
            (squeeze(sum(sum(sum(IC_probs_99,1),2),3)) + sum(curr_network(:)) ...
            - squeeze(sum(sum(sum(curr_IC_network_overlap,1),2),3)));
    end
    
    Group_QC.IC_in_network = IC_in_network;
    Group_QC.network_in_IC = network_in_IC;
    Group_QC.dice_IC_network = dice_IC_network;
    
    % Idea: select ICs that contribute toward networks:
    % For each network, sort ICs by dice . Loop through those calculating new
    % dice for each IC addition. Once new_dice < old_dice, stop. Output ICs
    % that contributed
    networks_ICs = cell(size(networks,4),1);
    networks_dice = zeros(size(networks,4),1);
    networks_IC_probs_99 = zeros(size(networks));
    networks_IC_probs_final = networks_IC_probs_99;
    for i1 = 1:size(networks,4)
        curr_dice = 0;
        curr_network_dices = dice_IC_network(:,i1);
        curr_network = networks(:,:,:,i1);
        [a, b] = sort(curr_network_dices, 'descend');
    
        % initialize
        curr_IC_prob_cum = zeros(size(networks(:,:,:,1)));
        go = 1; i2 = 1; network_ICs = [];
        while go > 0
            curr_IC_prob_cum = curr_IC_prob_cum + IC_probs_99(:,:,:,b(i2));
            curr_IC_and_network = logical(curr_IC_prob_cum .* curr_network); % A & B
            curr_IC_or_network = logical(curr_IC_prob_cum + curr_network); % A or B
            % dice = A&B / (AorB-A&B)
            dice_new = sum(curr_IC_and_network(:)) ./ sum(curr_IC_or_network(:));
            if dice_new <= curr_dice
                go = 0;
            else
                curr_dice = dice_new;
                network_ICs = [network_ICs, b(i2)]; % save a list of the ICs to use for each network
                updated_IC_prob_99 = int8(logical(curr_IC_prob_cum));
            end
            i2 = i2+1;
        end
        
        % writeout IC prob network file:
        networks_IC_probs_99(:,:,:,i1) = updated_IC_prob_99;
        networks_IC_probs_final(:,:,:,i1) = logical(updated_IC_prob_99 .* int8(curr_network)); % keep within network file
    
        networks_dice(i1) = curr_dice;
        networks_ICs{i1} = {network_ICs};
    end
    
    IC_probs_info.Datatype = 'single';
    IC_probs_info.ImageSize = [IC_probs_info.ImageSize(1:end-1), size(networks,4)];
    niftiwrite(cast(networks_IC_probs_99, 'single'), 'IC_probs_networks', IC_probs_info, 'Compressed', true)
    niftiwrite(cast(networks_IC_probs_final, 'single'), 'IC_probs_networks_final', IC_probs_info, 'Compressed', true)
    
    Group_QC.networks_dice = networks_dice;
    Group_QC.networks_ICs = networks_ICs;

end

save('group_qc.mat', 'Group_QC')
% When this is done, create your own table, with same image_numbers and
% image_names, with a list of what you are keeping versus not. Go through
% any you don't want to keep for qc reasons, and check if manual IC
% selection can be improved to help save it.
