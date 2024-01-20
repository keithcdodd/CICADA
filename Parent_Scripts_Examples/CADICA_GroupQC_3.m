% Run Group QC Analysis for a group of subjects
% 
% This can help you better evaluate how well CADICA did for each subject,
% and help make desicisions on data to keep or reject based on a number of
% parameters. It may also help you find if certain ICs of interest may have
% been missed persubject

% Run this on everyone first, and then you can use the excel documents to
% sort the qc values to find data to try manual IC adjustments for (if we
% are doing mainly auto CADICA) and to identify data that is simply
% bad/unsaveable. Put those in as the bad_data and the adjust prefixes

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%
%subj_ids = {'102', '107', '108'};
subj_ids = {'102', '103', '105', '106', '107', '108', '109', '112', '113', ...
    '114', '115', '116', '117', '118', '121', '122', '125', '126', '128', '129', '130', ...
    '134', '135', '138', '139', '140', '142', '143', '144', '145', '146', ...
    '147', '149', '153', '156', '158', '160', '161', '164', '165', '168', ...
    '169', '171', '172', '173', '174', '176', '178', '181', '184', '185', '186', '187'}; % List all subjects you want included as a group together
% subj_ids = {'102', '103'};
sess_ids = {'01', '02', '03', '04'}; % a single session id, might make sense for comparisons, but any you put in here will be included if they exist
task_ids = {'rest'}; % include runs if they exist, e.g., 'foodpics_run-01', but just do a single task name type - makes the most sense here for group comparisons
%bad_data_prefixes = {'sub-108_ses-02_task-rest', 'sub-144_ses-02_task-rest'}; % bad and not saveable from manual
bad_data_prefixes = {'sub-108_ses-02_task-rest', 'sub-144_ses-02_task-rest', 'sub-160_ses-01_task-rest', 'sub-172_ses-02_task-rest', 'sub-181_ses-02_task-rest'};
adjust_prefixes = {''}; % was seen as less than idea from auto, was adjusted manually
% adjust_prefixes = {'sub-102_ses-02_task_rest', 'sub-106_ses-02_task-rest', ...
%     'sub-108_ses-02_task-rest', 'sub-112_ses-01_task-rest', 'sub-112_ses-02_task-rest', ...
%     'sub-117_ses-03_task-rest', 'sub-117_ses-04_task-rest', 'sub-121_ses-02_task-rest', ...
%     'sub-122_ses-01_task-rest', 'sub-134_ses-01_task-rest', 'sub-134_ses-03_task-rest', ...
%     'sub-138_ses-02_task-rest', 'sub-138_ses-03_task-rest', 'sub-138_ses-04_task-rest', ...
%     'sub-139_ses-01_task-rest', 'sub-139_ses-02_task-rest', 'sub-139_ses-03_task-rest', ...
%     'sub-139_ses-04_task-rest', 'sub-142_ses-01_task-rest', 'sub-144_ses-02_task-rest', ...
%     'sub-147_ses-01_task-rest', 'sub-147_ses-02_task-rest', 'sub-153_ses-02_task-rest', ...
%     'sub-156_ses-04_task-rest', 'sub-158_ses-01_task-rest', 'sub-160_ses-01_task-rest', ...
%     'sub-160_ses-02_task-rest', 'sub-165_ses-01_task-rest', 'sub-165_ses-02_task-rest', ...
%     'sub-172_ses-02_task-rest', 'sub-173_ses-01_task-rest', 'sub-173_ses-02_task-rest', ...
%     'sub-174_ses-02_task-rest', 'sub-181_ses-01_task-rest', 'sub-181_ses-02_task-rest', ...
%     'sub-184_ses-02_task-rest', 'sub-185_ses-01_task-rest', 'sub-186_ses-02_task-rest'};
tr = 2;
redo_melodic = 0; % Whether or not to redo group melodic if it already exists - 0 is to not redo
gm_mni_prob_file = '/home/doddke/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_gm_tal_nlin_asym_09c_2mm.nii.gz'; % make sure same sampling as functionals
background_file = '/home/doddke/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_2mm.nii.gz'; % to use as background image in report for melodic
cadica_dir = '/data/images/awesome/data/bids_data/derivatives/cadica_2mm_stc_kd'; % 2mm resolution output space, slice time corrected, Keith Dodd
cadica_group_dir = '/data/images/awesome/data/bids_data/derivatives/cadica_2mm_stc_kd/groupqc_rest'; % folder you want to store the group QC information for this
network_file = ['/home/doddke/templates/network_template_', task_ids{1}, '.nii.gz'];
num_ICs = 20; % 20-30 is often a good number to break up general group ICs
% selecting which group of outputs you want to compare
CADICA_type = 'auto'; % manual or auto, manual is preferred if you have gone through it, need this to determine which folder you are using. Default is manual if not given.
smooth_filter_tag ='s_bp_'; % can be blank like '' or 's_' or 's_hp_' etc., if not blank, then it needs to end with underscore
reg_type = 'nonagg'; % either nonagg (nonaggressive - default, preferred) or agg (aggressive) 
%denoise_tag = ['CADICA_', CADICA_type, '_', reg_type]; % e.g., 'CADICA_manual_nonagg', 'CADICA_auto_agg', 'orig', '8p', '9p', etc., whatever text is between the smoothing_filtering_ and _bold.nii.gz
denoise_tag = ['8p'];
% denoise_type_tag needs to be an exact match to your files. This allows
% you to make your own type of denoising too, but naming needs to follow
% similar conventions
% compare tags is to decide what folders of comparison images you want to
% group together. May be suggested to do at least 8p and
% CADICA_auto_nonagg, assuming you are running the manual nonaggressive as
% your main data.
%compare_tags = {'8p', '9p', 'CADICA_auto_nonagg'}; % find in your qc folder the .mat files, the tag you are looking for is inbetween 'vs_' and '_qc'. '8p' or '9p' would be great standard comparisons. So would an auto CADICA
compare_tags = {'8p', '9p'}; % if doing auto
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For FSL to run all happy, we need to set the LD_LIBRARY_PATH Differently.
setenv('LD_LIBRARY_PATH', '/usr/lib64/lib64/libstdc++.so.6:/usr/local/MATLAB/R2022b/bin/matlab')
addpath('/home/doddke/code/CADICA')

% OK, we need to go through each file of interest (per subject, session,
% task) and then load the 3D Signal map, alongside save number of ICs kept,
% percent ICs kept, and new folder for all this and to store the manual vs
% aggressive (or which ever thing we care about the most) for quick
% comparisons as well.

% convert gm probability into a mask at .67
gm_mni_prob = niftiread(gm_mni_prob_file);
gm_mni_thresh = gm_mni_prob_file > 0.67;
gm_mni_thresh_info = niftiinfo(gm_mni_prob_file);
gm_mni_thresh_info.Datatype = 'uint8'; % in case I want to write it out as a file later


% find the correct CADICA cleaned folder tag, default is manual
if strcmp(CADICA_type, 'auto') ~= 1
    CADICA_type = 'manual';
end

CADICA_type_default = CADICA_type;

% Create cadica_group_dir, if it does not already exist:
if not(isfolder(cadica_group_dir))
    mkdir(cadica_group_dir)
end

% Create a folder that is more specific within there
qc_specific_fol = [cadica_group_dir, '/', smooth_filter_tag, denoise_tag];
if not(isfolder(qc_specific_fol))
    mkdir(qc_specific_fol)
end

data_fol = [qc_specific_fol, '/data'];
if not(isfolder(data_fol))
    mkdir(data_fol)
else
   rmdir(data_fol, 's') % remove data folder and contents if it already exists
end

% Make separate photo folders, so it is easier to just scroll through them
% all in the future
for h = 1:length(compare_tags)
    % create folders for comparison photos to zip through too
    photo_fol = [qc_specific_fol, '/qc_photos_', compare_tags{h}];
    if not(isfolder(photo_fol))
        mkdir(photo_fol)
    end
end

% don't forget to initialize with zeros a 4D signaltonoise file to save on
% time
m = 1; % counting purposes, nothing major or wild

for j = 1:length(subj_ids)

    fprintf('\n')
    fprintf(['Running for subject ', subj_ids{j}, '\n'])
    % check for folder
    curr_subj_fol = [cadica_dir, '/sub-', subj_ids{j}];
    if ~isfolder(curr_subj_fol)
        fprintf(['No subject folder found for ', subj_ids{j}, '. Moving on to next one...\n'])
        continue
    end

    for k = 1:length(sess_ids)
        
        fprintf(['Running for session ', sess_ids{k}, '\n'])
        % check for folder
        curr_sess_fol = [curr_subj_fol, '/ses-', sess_ids{k}];
        if ~isfolder(curr_sess_fol)
            fprintf(['No sess folder found for ', sess_ids{k}, '. Moving on to next one...\n'])
            continue
        end

        for l = 1:length(task_ids)
            fprintf(['Running for task ', task_ids{l}, '\n'])
            curr_task_fol = [curr_sess_fol, '/', task_ids{l}];
        if ~isfolder(curr_task_fol)
            fprintf(['No task folder found for ', task_ids{l}, '. Moving on to next one...\n'])
            continue
        end

            cd(curr_task_fol)

            % OK, now grab the things you need/want:
            % The actually cleaned file, in cleaned_manual folder:
            prefix = ['sub-', subj_ids{j}, '_ses-', sess_ids{k}, '_task-', task_ids{l}];
            suffix = 'bold.nii.gz';
            % check if this is a subject where you want to use the manually
            % adjusted, or not include it
            % To Do: figure out what you want to do here then
            adjusted = 0; % 1 if we adjusted it
            bad_data = 0; % if it was bad data
            CADICA_type = CADICA_type_default;
            denoise_tag = ['CADICA_', CADICA_type, '_', reg_type];
            if sum(contains(adjust_prefixes, prefix)) > 0
                CADICA_type = 'manual';
                denoise_tag = ['CADICA_', CADICA_type, '_', reg_type];
                adjusted = 1;
                bad_data = 0;
            elseif sum(contains(bad_data_prefixes, prefix)) > 0
                bad_data = 1; % we can still run it with the manual version, why not
                CADICA_type = 'manual';
                denoise_tag = ['CADICA_', CADICA_type, '_', reg_type];
                adjusted = 1;
            end

             % check for cleaned folder
            cleaned_folname = ['cleaned_' CADICA_type];
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
            
            selection_fol = [curr_task_fol, '/ic_', CADICA_type, '_selection'];
            curr_signalandnoise_file = [selection_fol, '/', curr_signalandnoise_overlap_filename];

            curr_signalandnoise_mask = niftiread(curr_signalandnoise_file);

            % Calculate overlap with GM mask
            if m == 1
                flirt_command = ['flirt -ref ', curr_signalandnoise_file, ' -in ', gm_mni_prob_file, ' -out ', qc_specific_fol, '/gm_mni_prob_resampled.nii.gz -usesqform -applyxfm'];
                [status, cmdout_flirt] = system(flirt_command, '-echo');
                gm_mni_prob_resampled = niftiread([qc_specific_fol, '/gm_mni_prob_resampled.nii.gz']);
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
            percent_ICs_kept = num_ICs_kept / num_ICs_total;

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
            DOF_estimate_if_standard_bp = percent_variance_kept * percent_freq_kept_if_standardbp * numvolumes;

            qc_group_array(:,m) = {m, curr_filename, subj_ids{j}, sess_ids{k}, task_ids{l}, adjusted, bad_data, gm_signal_dice, num_ICs_kept, ...
                num_ICs_total, percent_ICs_kept, percent_variance_kept, ...
                percent_freq_kept, numvolumes, DOF_estimate_final, DOF_estimate_if_standard_bp};

            data_files{m} = curr_copied_denoised_file;

            % finally copy over the photos to their respective folders
            for h = 1:length(compare_tags)
                photo_fol = [qc_specific_fol, '/qc_photos_', compare_tags{h}];
                photo_name = [prefix, '_', denoise_tag, '_vs_', compare_tags{h}, '_qc_plots.png'];
                photo_file = [curr_task_fol, '/qc/', photo_name];
                copyfile(photo_file, photo_fol)
            end
            m = m+1; % increment m

        end
    end
end

% go to general group folder for this specific set
cd(qc_specific_fol)

% collapse funcmask into 3D with mean
mean_funcmask = mean(funcmask, 4); 
funcmask_info = niftiinfo(curr_mask_file); % will be 3D instead of 4D, will need 
niftiwrite(cast(mean_funcmask, 'uint8'), 'mean_funcmask_thresh', funcmask_info, 'Compressed', true) % 0-1, it should be the percent of funcmasks that overlapped with voxel
funcmask_info.Datatype = 'single';
niftiwrite(cast(mean_funcmask, 'single'), 'mean_funcmask_prob', funcmask_info, 'Compressed', true) % 0-1, it should be the percent of funcmasks that overlapped with voxel

% Create 4d funcmask too so it is easy to search through
funcmask_info.Datatype = 'uint8';
funcmask_info.ImageSize = [funcmask_info.ImageSize, size(signalandnoise_mask,4)];
funcmask_info.PixelDimensions = [funcmask_info.PixelDimensions, 1];
niftiwrite(cast(funcmask, 'uint8'), 'funcmask_Group', funcmask_info, 'Compressed', true)

% OK, now load and create 4D signaltonoise file and header, and write to
% the group directory
signal_prob_info = niftiinfo(curr_signalandnoise_file); % will be 3D instead of 4D, will need to edit to 4D
signal_prob_info.ImageSize = [signal_prob_info.ImageSize, size(signalandnoise_mask,4)];
signal_prob_info.PixelDimensions = [signal_prob_info.PixelDimensions, 1]; % idk if 1 is "right", but it makes some sense/should work! 4 dimension is just each image

niftiwrite(cast(signalandnoise_mask, 'single'), 'SignalandNoise_Group', signal_prob_info, 'Compressed', true) % 1 is high signal, -1 is high noise-labeled only
%niftiwrite(cast(signalandnoise_mask, 'single'), 'Signal_Group', signal_prob_info, 'Compressed', true) % 1 is high signal, -1 is high noise-labeled only


% cell to table
qc_table_labels = {'image_number', 'image_names', 'subject', 'session', 'task', 'manually_adjusted', 'removed_bad'...
    'gm_signal_dice', 'number_kept_ics', 'number_total_ics', ...
                'fraction_kept_ics', 'fraction_signal_variance_kept', 'fraction_frequency_kept', 'numvolumes', ...
                'dof_estimate_final', 'dof_estimate_final_if_standard_resting_bp'};
final_qc_table = cell2table(qc_group_array', 'VariableNames', qc_table_labels);
writetable(final_qc_table, 'group_qc_table.csv') % can use this and sort by gm_signal and signal_gm proportions and also number_kept_ics to help find the bad or adjustable data

save('qc_table_vals.mat', 'final_qc_table')

% Run Group MELODIC
% First, make text file of the image names
% will need a text file where separated by commas (no spaces) and then in
% the melodic call $(cat )
% Make sure to not include data files that were designated as bad enough to
% remove
writecell(data_files(~final_qc_table.removed_bad), 'image_names.txt')

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
        [status, cmdout_ICA] = system(Group_ICA_command, '-echo');
        fprintf('Done with Group ICA (Melodic)\n\n')
    end
elseif redo_melodic == 1
    fprintf(['Running: ', Group_ICA_command, '\n'])
    [status, cmdout_ICA] = system(Group_ICA_command, '-echo');
    fprintf('Done with Group ICA (Melodic)\n\n')
end

% create ICprobabilities combined in melodic folder (good for making masks
% later at 100% probability, for example).
probmapnames = '"$(ls ./melodic/stats/probmap_* | sort -V)"';
command = ['fslmerge -t ./melodic/ICprobabilities.nii.gz ', probmapnames];
[status, cmdout] = system(command, '-echo');
% And load this now, 4th dimension is the ICs
IC_probs = niftiread('/melodic/ICprobabilities.nii.gz');
IC_probs_info = niftiinfo('/melodic/ICprobabilities.nii.gz');
IC_probs_95 = logical(IC_probs > 0.949); % keeps only best overlap

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
% because matlab does not love a 4th dimension, the code might have to get
% weird...
% Will want to save values into some sort of chart per IC and per network
for i1 = 1:size(networks,4)
    curr_network = networks(:,:,:,i1);
    curr_IC_network_overlap = logical(IC_probs_95 == curr_network); % Try to mask in 4th dimension
    IC_in_network = curr_IC_network_overlap ./ sum(IC_probs_95,4);
    network_in_IC = curr_IC_network_overlap ./ sum(curr_network(:));
    dice_IC_network = curr_IC_network_overlap ./ (sum(IC_probs_95,4) + sum(curr_network(:)) - curr_IC_network_overlap);
end

% When this is done, create your own table, with same image_numbers and
% image_names, with a list of what you are keeping versus not. Go through
% any you don't want to keep for qc reasons, and check if manual IC
% selection can be improved to help save it.
