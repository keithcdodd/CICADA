% script to run group qc of cicada files
% This will all be based on the same .csv you used for initial cicada runs.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%
file_tag = 'auto'; % a tag that is unique to the type of denoised file. E.g., 'auto', 'manual', '8p', '9p', etc.
cicada_csv = '/Users/keithdodd/Work/awesome/documents/cicada_runs/cicada_runs.csv'; % As described above, change as needed
rerun_melodic_regardless = 0; % 0 will not rerun melodic if it is found in default space. 1 will rerun regardless
cicada_wrapper_dir = '/Users/keithdodd/Work/code/CICADA/wrappers';
task_name = 'rest'; % should only be on one task for group qc
redo_melodic = 1; % 0 is to not redo if it is already done. 1 redoes it even if already done.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% need to make sure CICADA folder and subfolders are added to path!
CICADA_group_scripts_dir = fileparts(mfilename('fullpath')); % this gives current script path
addpath(CICADA_group_scripts_dir); % current location to path 

% make sure everything is read as a cell array of char vector from the .csv
% for consistency
opts = detectImportOptions(cicada_csv);
opts = setvartype(opts, 'char');  

cicada_runs_table = readtable(cicada_csv, opts); % now read in the actual inputs!
cicada_runs = cicada_runs_table(strcmp(cicada_runs_table.task_name, task_name),:);
num_runs = size(cicada_runs,1);

fprintf('###########################################\n')
fprintf(['Running for ', num2str(num_runs), ' number of runs!\n'])

% read variables, make sure excel has these columns
cicada_dirs = cicada_runs.cicada_dir;
output_dirs = cicada_runs.group_qc_dir; % Add this to your excel
sub_ids = cicada_runs.sub_id;
ses_ids = cicada_runs.ses_id;
task_names = cicada_runs.task_name;
initial_outliers = cicada_runs.initial_outlier; % data that was unusable from the get-go
denoising_outliers = cicada_runs.denoising_outlier; % data where denoising was unable to save the data, or is deemed to throw out for other reasons
% NOTE: An example for using denoising_outlier might be when the data
% exceeds movement thresholds, or if CICADA marked outliers, and you want
% to apply this same outlier list (manually) here. The point is, you, or a
% program, determines outliers due to thresholds and this is where you can
% manually mark that. The data will still be a part of group_qc, but for
% future analyses, you can use this list to not include them.
manually_adjusteds = cicada_runs.manually_adjusted; % data you have tried to save with manual ICA selection
task_event_files = cicada_runs.task_events_file;

% check that group_qc_dir and cicada_dir are the same for everyone (it should be)
if ~(sum(strcmp(output_dirs, output_dirs{1})) == length(output_dirs))
    fprintf('Group QC directory is not the same for everyone!\n')
    return;
end

if ~(sum(strcmp(cicada_dirs, cicada_dirs{1})) == length(cicada_dirs))
    fprintf('CICADA directory is not the same for everyone!\n')
    return;
end

% set up folders for outputs:
% Create output_dir, if it does not already exist:
output_dir = [output_dirs{1}, '/', task_name]; % because they should ALL be the same, and specified to task
if not(isfolder(output_dir))
    mkdir(output_dir)
end

% create a data folder, remove old one if it exists
data_fol = [output_dir, '/data'];
if not(isfolder(data_fol))
    mkdir(data_fol)
else
   rmdir(data_fol, 's') % remove data folder and contents if it already exists
   mkdir(data_fol) % and now make new folder
end


% Make separate photo folders, so it is easier to just scroll through them
% all in the future. There is only 8p compare if auto, but 8p and auto if
% manual 

% go into first participant's info to extract the correct naming:
first_task_dir = [cicada_dirs{1}, '/sub-', sub_ids{1}, '/ses-', ses_ids{1}, '/', task_name];
compare_image_info = dir([first_task_dir, '/qc/sub*ses*task*', file_tag, '*vs*_qc_plots.jpg']); % specific to file tag of interest, can be multiples

% make sure that the dir call actually found the file(s) of images that
% compare to your cleaned file of interest (given the file tag)
if size(compare_image_info,1) == 0
    fprintf(['Could not find a compare file match at task directory at ', first_task_dir, '/qc/ \nIs your file_tag correct, or does the task dir not exist?\n'])
    return;
end

% Ok, NOW this should make appropriate photo folders
for h = 1:size(compare_image_info,1)
    compare_tag = extractBetween(compare_image_info(h).name, [task_name, '_'], '_qc_plots.jpg');
    compare_tag = compare_tag{:};
    photo_fol = [output_dir, '/qc_photos_', compare_tag];
    qc_photo_fols{h} = photo_fol; % for use for later for copy and past
    if not(isfolder(photo_fol))
        mkdir(photo_fol)
    else
        rmdir(photo_fol, 's')
        mkdir(photo_fol)
    end
end

% Now go into first cleaned_dir in order to extract if these are CICADA
% files or not
first_cleaned_dir = [first_task_dir, '/cleaned'];
first_image_info = dir([first_cleaned_dir, '/*', file_tag, '*.nii.gz']); % specific to file tag of interest, should be only ONE file

if size(first_image_info,1) ~= 1
    fprintf(['File tag', file_tag, ', is either not specific enough to only grab one file, or grabs no files. \n'])
    return;
end

% grab tr, which can be relevant for melodic 
first_image_nifi_info = niftiinfo([first_image_info.folder, '/', first_image_info.name]);
tr = first_image_nifi_info.PixelDimensions(4); % for use later with melodic

if contains(first_image_info.name, 'CICADA')
    cicada = 1; % mark if these will be CICADA files or not, important for qc tables
else
    cicada = 0;
end

% the following block is needed for later (like making group melodic look
% easy to follow with a nice mni background image):
cicada_group_qc_dir = fileparts(mfilename('fullpath')); % this gives current script path
cd([cicada_group_qc_dir, '/..'])
cicada_path = pwd;
gm_mni_prob_file = [cicada_path, '/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_gm_tal_nlin_asym_09c.nii.gz'];
background_file = [cicada_path, '/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c.nii.gz'];
network_file = [cicada_dirs{1}, '/templates/network_template_', task_name, '.nii.gz'];

% convert gm probability into a mask at .67
gm_mni_prob = niftiread(gm_mni_prob_file);
gm_mni_thresh = gm_mni_prob_file > 0.67;
gm_mni_thresh_info = niftiinfo(gm_mni_prob_file);
gm_mni_thresh_info.Datatype = 'uint8'; % in case I want to write it out as a file later

% initialize struct to store all information:
Group_QC = struct;

% because not all instances might work, we should keep our own counter
% within the for loop 
m = 1;
bd = 1;

for idx = 1:num_runs
    cicada_dir = cicada_dirs{idx}; % should be the same for all
    output_dir = output_dirs{idx}; % should be the same for all
    sub_id = sub_ids{idx};
    ses_id = ses_ids{idx};
    task_name = task_names{idx}; % should be the same for all 
    bad_data = initial_outliers{idx}; % needs to be '1' for bad data, or '0' for OK data
    manually_adjusted = manually_adjusteds{idx}; % needs to be '1' for manually adjsuted, '0' for not
    denoising_outlier = denoising_outliers{idx};
    task_event_file = task_event_files{idx};

    fprintf(['\nRunning sub ', sub_id, ' ses ', ses_id, ' task ', task_name, '\n'])

    % if bad data, whole run needs to be skipped and not included in qc
    % (use continue)
    if ~strcmp(bad_data, '0')
        % this is bad data from the get go (recognized previously, like
        % with mriqc, or at time of scan), do not include in group qc
        % analysis
        bad_data_prefixes{bd} = ['sub-', sub_id, '_ses-', ses_id, '_task-', task_name];
        fprintf(['   Sub ', sub_id, ' ses ', ses_id, ' task ', task_name, ' was marked as bad data. Skipping...\n'])
        continue;
    end

    % Make sure task directory exists
    task_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];
    if ~isfolder(task_dir)
        fprintf('   Cannot find task directory at ', task_dir, '. Skipping...\n')
        continue;
    end

    % cicada type determines if we default from pulling from cleaned_auto
    % or cleaned_manual. Manually adjusted ('1') would mean pull from manual too.
    cleaned_dir = [task_dir, '/cleaned'];

    % then, test if desired cleaned dir exists
    if ~isfolder(cleaned_dir)
        fprintf(['   Cannot find designated cleaned_dir at ', cleaned_dir, '. Skipping...\n'])
        continue;
    end

    % if manually adjusted, will need to pull from manual for this task
    % instance

    % first grab general file tag and see if it contains CICADA in it:
    test_cleaned_file_info = dir([cleaned_dir, '/*', file_tag, '*.nii.gz']);

    % make sure only one file is grabbed
    if size(test_cleaned_file_info,1) ~= 1
        fprintf('Either your file_tag is not specific enough (grabs more than one file) OR no files match with file_tag\n')
        return;
    end

    if contains(test_cleaned_file_info.name, 'CICADA')
        % This confirms we are working with CICADA data, and thus should
        % check if manually adjusted is the preferred usage
        if ~strcmp(manually_adjusted, '0')
            % this was manually adjusted, grab manual version
            cleaned_file_info = dir([cleaned_dir, '*CICADA*manual*.nii.gz']);

            % make sure only one file is grabbed
            if size(cleaned_file_info,1) ~= 1
                fprintf('Do you have more than one CICADA*manual*.nii.gz file somehow?\n')
                return;
            end
        else
            % else, use file_tag
            cleaned_file_info = test_cleaned_file_info;
        end
    else
        % this is not a cicada file, so we do not need to worry about
        % manual CICADA adjustments
        cleaned_file_info = test_cleaned_file_info;
    end

    cleaned_file = [cleaned_file_info.folder, '/', cleaned_file_info.name];

    
    % OK, now FINALLY loop through cleaned_dir files, and run qc for all of them!
    % then finally individually grab relevant qc
    [cleaned_data, data_mask, signalandnoise_overlap, qc_table, qc_corrs_table, qc_photo_paths] = cicada_get_qc(cleaned_file);

    % then we need to be adding this information into the group qc stuff
    % (Group_QC struct)
    cleaned_data_list{m} = cleaned_data;
    data_mask_list{m} = data_mask;
    signalandnoise_overlap_list{m} = signalandnoise_overlap;
    denoising_outlier_list{m} = strcmp(denoising_outlier, '1');
    task_event_file_exist_list{m} = ~isempty(task_event_file);

    if m == 1
         % if it is first instance, then group qc table is just subect qc
         % table
         group_qc_table = qc_table;
         group_qc_corrs_table = qc_corrs_table;
    else
        % otherwise, tack it on to the end
        group_qc_table = [group_qc_table; qc_table];
        group_qc_corrs_table = [group_qc_corrs_table; qc_corrs_table];
    end

    % And then move the qc photos to their respective folders!
    compare_image_info = dir([task_dir, 'sub*ses*task*', file_tag, '*vs*_qc_plots.jpg']); % specific to file tag of interest, can be multiples
    for h = 1:size(compare_image_info,1)
        curr_compare_image_file = [compare_image_info(h).folder, '/', compare_image_info(h).name]; % use later for copy and paste
        compare_tag = extractBetween(compare_image_info(h).name, [task_name, '_'], '_qc_plots.jpg');
        compare_tag = compare_tag{:};

        photo_fol = [output_dir, '/qc_photos_', compare_tag];
        copyfile(curr_compare_image_file, photo_fol)
    end
   
    m = m+1; % increment successful counter
end

% record the bad data that was not even looked at, for easy reference later:
Group_QC.bad_data_prefixes = bad_data_prefixes;

% go to general group folder for this specific set
cd(output_dir)

% Calculate outliers based on the final qc table. Do it in three ways and
% add it to the qc table (only if this is CICADA though!):
final_qc_table = group_qc_table;

% Now we can calculate outliers, but some are only calculated in this
% manner IF it is CICADA data
if cicada == 1
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
    if sum(cell2mat(task_event_file_exist_list)) > 0 % so there are task event files!
        fprintf('Found task event files, so using best task power overlap.\n')
        Group_QC.low_besttask_power_overlap = (isoutlier(final_qc_table.besttask_power_overlap, "median")) & ...
        (final_qc_table.besttask_power_overlap < mean(final_qc_table.besttask_power_overlap));
        Group_QC.low_power_overlap = Group_QC.low_general_power_overlap & Group_QC.low_besttask_power_overlap;
    else
        %else it is resting state
        fprintf('No Task event files. Assumed to be resting state.\n')
        Group_QC.low_power_overlap = Group_QC.low_general_power_overlap;
    end
    
    
    % DVARS_Corr
    Group_QC.high_DVARS_corr = (isoutlier(abs(final_qc_table.DVARS_GM_median), "median")) & (abs(final_qc_table.DVARS_GM_median) > median(abs(final_qc_table.DVARS_GM_median)));
    
    % FD_Corr
    Group_QC.high_FD_corr = (isoutlier(abs(final_qc_table.FD_GM_median), "median")) & (abs(final_qc_table.FD_GM_median) > median(abs(final_qc_table.FD_GM_median)));
    
    
    % combine to get overall outliers
    Group_QC.cicada_outliers = logical(Group_QC.low_number_total_ics + ...
        Group_QC.low_fraction_signal_variance_kept + Group_QC.low_Signal + ...
        Group_QC.low_Smoothing + Group_QC.low_power_overlap + Group_QC.low_gm_coverage + Group_QC.low_gm_dice);

    % add to final qc table
    final_qc_table.low_number_total_ics = Group_QC.low_number_total_ics;
    final_qc_table.low_fraction_signal_variance_kept = Group_QC.low_fraction_signal_variance_kept;
    final_qc_table.low_Signal = Group_QC.low_Signal;
    final_qc_table.low_Smoothing = Group_QC.low_Smoothing;
    final_qc_table.low_power_overlap = Group_QC.low_power_overlap;
    final_qc_table.low_gm_coverage = Group_QC.low_gm_coverage;
    final_qc_table.high_DVARS_corr = Group_QC.high_DVARS_corr;
    final_qc_table.high_FD_corr = Group_QC.high_FD_corr;
    final_qc_table.cicada_outliers = Group_QC.cicada_outliers;   
end

% Now, get liberal outliers (not restrictive)
Group_QC.liberal_outliers = final_qc_table.mean_FD > 0.55;

% Finally, conservative outliers: median FD
Group_QC.conservative_outliers = (final_qc_table.mean_FD > 0.25) | ...
    (final_qc_table.Percent_FD_gt_point2mm > 20) | ...
    (final_qc_table.AnyFD_gt_5mm);

% And record Denoising outlier (which was manually selected by the user)
Group_QC.denoising_outliers = cell2mat(denoising_outlier_list)';

% add to final qc table
final_qc_table.liberal_outliers = Group_QC.liberal_outliers;
final_qc_table.conservative_outliers = Group_QC.conservative_outliers;
final_qc_table.denoising_outliers = Group_QC.denoising_outliers;

writetable(final_qc_table, 'group_qc_table.csv') % can use this and sort by gm_signal and signal_gm proportions and also number_kept_ics to help find the bad or adjustable data
Group_QC.final_qc_table = final_qc_table;

% Now work on the images:

% remove outliers from the rest of analysis
if cicada == 1
    image_keep = ~Group_QC.cicada_outliers & ~Group_QC.denoising_outliers;
else
    image_keep = ~Group_QC.denoising_outliers;
end

% save image_keep so it is easy in the future to know what images were left
% out
Group_QC.Images_Used = image_keep;

% merge signalandnoise_overlap into a 4D file, if CICADA
if cicada == 1
    merge_signalandnoise_overlap_command = ['fslmerge -a ', output_dir, '/signal_noise_overlaps.nii.gz ', strjoin(signalandnoise_overlap_list(image_keep))];
    [~, ~] = call_fsl(merge_signalandnoise_overlap_command);
end

% merge func masks as well!
merge_funcmasks_command = ['fslmerge -a ', output_dir, '/funcmasks.nii.gz ', strjoin(data_mask_list(image_keep))];
[~, ~] = call_fsl(merge_funcmasks_command);

% and calculate a Tmin of data masks (where there is funcmask overlap
% across all subects)
funcmask_overlap_command = ['fslmaths ', output_dir, '/funcmasks.nii.gz -Tmin ', output_dir, '/group_funcmask.nii.gz'];
[~, ~] = call_fsl(funcmask_overlap_command);


% Run Group MELODIC
% First, make text file of the image names
% will need a text file where separated by commas (no spaces) and then in
% the melodic call $(cat )
% Make sure to not include data files that were designated as bad enough to
% remove
data_files = cleaned_data_list;
all_data_files = data_files;
data_files = all_data_files(image_keep);
writecell(data_files, 'image_names.txt')

% Resample mni to be the background image for melodic
flirt_command = ['flirt -ref group_funcmask.nii.gz -in ', background_file, ' -out mni_resampled.nii.gz -usesqform -applyxfm'];
[~, ~] = call_fsl(flirt_command);

% If you want to specify how many ICs to make:
%Group_ICA_command = ['melodic --in=$(cat image_names.txt) --outdir=melodic --Ostats -m mean_funcmask_thresh.nii.gz --bgimage=mni_resampled.nii.gz -d ', num2str(num_ICs), ' --nobet --mmthresh=0.5 --report --tr=', num2str(tr)];
% Vs.: let Melodic decide on number of ICs (preferred)
Group_ICA_command = ['melodic --in=$(cat image_names.txt) --outdir=melodic --Ostats -m group_funcmask.nii.gz --bgimage=mni_resampled.nii.gz --nobet --mmthresh=0.5 --report --tr=', num2str(tr)];

if redo_melodic == 0
    if not(isfolder('melodic'))
        fprintf(['Running: ', Group_ICA_command, '\n'])
        fprintf('Done with Group ICA (Melodic)\n\n')
    end
elseif redo_melodic == 1
    fprintf(['Running: ', Group_ICA_command, '\n'])
    [~, ~] = call_fsl(Group_ICA_command);
    fprintf('Done with Group ICA (Melodic)\n\n')

    % create ICprobabilities combined in melodic folder (good for making masks
    % later at 100% probability, for example).
    % first get to qc folder in bash
    out = system(['cd ', output_dir], '-echo');
    [~, str] = system('ls ./melodic/stats/probmap_* | sort -V');
    probmapnames = regexprep(str, '\s+', ' '); % convert newlines to spaces
    command = ['fslmerge -t ./melodic/ICprobabilities.nii.gz ', probmapnames];
    [~, ~] = call_fsl(command);

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
