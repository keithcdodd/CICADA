function cicada_group_qc(cicada_home, group_qc_home, task_name, output_dirname, file_tag, smoothing_kernel, fpass, detrended_degree, redo_melodic, sub_ids, ses_ids, excludes, outliers, adjusteds, task_event_files)
% function to run group qc
% cicada_home: is the cicada home directory, the general home input folder
% group_qc_home: the general group qc output dir
% output_dirname: what name do you want for your output directory? This
% will end up like the following: group_qc_home/task_name/output_dirname
% task_name: The task id for the group qc data, e.g. 'visual_run-01', or
% 'rest'
% sub_ids, ses_ids, excludes, outliers, adjusteds, and
% task_event_files should all be cell arrays of characters of the same number of
% rows. This function will loop through each row of them all and match them
% for the call. 
% file_tag: needs to be unique to the type of denoised/cleaned file you are
% putting together in group qc. e.g., '_auto_', '_manual_', '_8p_'

% smoothing_kernel: gaussian smoothing kernel FWHMx size. Default is same
% as current size of functional (can do default by giving value of -1). QC
% plots are most accurate if no smoothing (or minimal smoothing) is
% applied. Smoothing will otherwise recorrelate nearby voxels and noise
% profiles in QC plots will also bleed more together.

% fpass is Hz for bandpassing. If used, 0.008 to 0.15 is recommended.
% Default is no bandpassing []. Lower and higher frequencies should already
% be diminished by ICA denoising. 

% detrended degree is the degree of polynomial to detrend at. Detrending to
% the second polynomial is the default if this is input incorrectly. 

% redo_melodic: whether or not to rerun the group melodic if it was already
% run. 0 to not rerun, 1 to rerun regardless

% sub_ids e.g., {'102', '102', '103'}
% ses_ids e.g., {'01', '02', '01'}
% excludes (either 0 or 1 for exclude) e.g., {'0', '1', '0'}
% outliers (either 0 or 1 for outlier) e.g., {'0', '1', '0'}
% adjusteds (either 0 or 1 for manually IC adjusted) e.g., {'0', '1', '0'}
% task_event_files should contain paths to the task_event_file(s)
% (optional), also as a cell array

% if no task_event_file_list is given (e.g., it is all resting state) then
% create an empty cell array of chars the same length as sub_id_list
if ~exist('task_event_files', 'var') || isempty(task_event_files)
    task_event_files = cell(size(sub_ids));
    task_event_files(:) = {''}; 
end

% now check that all lists are the same length, and that all are lists of
% char arrays
if ~isequal(length(sub_ids), length(ses_ids), length(excludes), length(outliers), length(adjusteds), length(task_event_files))
    fprintf('Your lists (sub_ids, ses_ids, excludes, outliers, adjusteds, task_event_files) exist at different lengths... quitting\n')
    return
end

test_class = {''};

if ~isequal(class(test_class), class(sub_ids), class(ses_ids), class(excludes), class(outliers), class(adjusteds), class(task_event_files))
    fprintf('Your lists are not all cell arrays \n')
    return
end

if ~isequal(class(test_class{1}), class(sub_ids{1}), class(ses_ids{1}), class(excludes{1}), class(outliers{1}), class(adjusteds{1}), class(task_event_files{1}))
    fprintf('Your lists do not contain all chars (should be cell arrays of characters!) \n')
    return
end

% if no smoothing kernel exists, make it 3mm (small, but helpful)
if ~exist('smoothing_kernel', 'var') || ~isnumeric(smoothing_kernel) || smoothing_kernel < 0
    fprintf('Default 3mm smoothing gaussian kernel will be applied\n')
    smoothing_kernel = 3;
end

% if no detrended degree exists, make it 2
if ~exist('detrended_degree', 'var') || ~isnumeric(detrended_degree)
    fprintf('Default 2nd degree detrending will be applied\n')
    detrended_degree = 2;
end

% set up folders for outputs:
% Create output_dir, if it does not already exist:
output_dir = [group_qc_home, '/', task_name, '/', output_dirname]; % because they should ALL be the same, and specified to task
if not(isfolder(output_dir))
    mkdir(output_dir)
end

if (~exist('fpass', 'var') == 1) || (~isa(fpass, 'double') == 1) || isempty(fpass)
    fpass = []; % if bandpass frequencies are not readable, assume no bandpass applied
end

% Make separate photo folders, so it is easier to just scroll through them
% all in the future. There is only 8p compare if auto, but 8p and auto if
% manual 

% go into first participant's info to extract the correct naming, just compare to 8p here for ease:
first_task_dir = [cicada_home, '/sub-', sub_ids{1}, '/ses-', ses_ids{1}, '/', task_name];
compare_image_info = dir([first_task_dir, '/qc/sub*ses*task*', file_tag, '*vs*8p*qc_plots.jpg']); % specific to file tag of interest, can be multiples

comparison_only = 0; % to change in following if else statement
orig_only = 0; % to change in following if else statement
% make sure that the dir call actually found the file(s) of images that
% compare to your cleaned file of interest (given the file tag)
if size(compare_image_info,1) == 0
    % check if it exists as a compare only instead
    compare_image_compareonly_info = dir([first_task_dir, '/qc/sub*ses*task*vs*', file_tag, '*qc_plots.jpg']); 
    
    if size(compare_image_compareonly_info,1) == 0
        % check if it exists as orig only
        if matches(file_tag, 'orig')
            orig_only = 1; % so now we will not use a photo directory for orig only
        else
            fprintf(['Could not find a compare file match at task directory at ', first_task_dir, '/qc/ \nIs your file_tag correct, or does the task dir not exist?\n'])
            return;
        end
    else
        % OK, so the group qc is being run on a file type that has only
        % ever been used as a comparison
        comparison_only = 1;
    end
end

% Ok, NOW this should make appropriate photo folders. If it is comparison
% only or orig only, then do not need (e.g., 8p)
if comparison_only == 0 && orig_only == 0
    for h = 1:size(compare_image_info,1)
    compare_tag = extractBetween(compare_image_info(h).name, [task_name, '_'], '_qc_plots.jpg');
    compare_tag = compare_tag{:};
    photo_fol = [output_dir, '/qc_photos_', compare_tag];
    qc_photo_fols{h} = photo_fol; % for potential use for later for copy and paste
        if not(isfolder(photo_fol))
            mkdir(photo_fol)
        else
            rmdir(photo_fol, 's')
            mkdir(photo_fol)
        end
    end
end

% Create a data folder, don't delete old one
data_dir = [output_dir, '/data'];
if isfolder(data_dir)
    rmdir(data_dir, 's')
end
mkdir(data_dir)

% Now go into first cleaned_dir in order to extract if these are CICADA
% files or not
first_cleaned_dir = [first_task_dir, '/cleaned'];

% a fail safe for funcmasks that could be lurking inside
if ~isempty(dir([first_cleaned_dir, '/funcmask*']))
    movefile([first_cleaned_dir, '/funcmask*'], [first_cleaned_dir, '/..'])
end

first_image_info = dir([first_cleaned_dir, '/*', file_tag, '*.nii.gz']); % specific to file tag of interest, should be only ONE file

if size(first_image_info,1) ~= 1
    fprintf(['File tag ', file_tag, ', is either not specific enough to only grab one file, or grabs no files. \n'])
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
network_file = [cicada_home, '/templates/network_template_', task_name, '.nii.gz'];

% convert gm probability into a mask at .67
%gm_mni_prob = niftiread(gm_mni_prob_file);
%gm_mni_thresh = gm_mni_prob_file > 0.67;
gm_mni_thresh_info = niftiinfo(gm_mni_prob_file);
gm_mni_thresh_info.Datatype = 'uint8'; % in case I want to write it out as a file later

% initialize struct to store all information:
Group_QC = struct;

% because not all instances might work, we should keep our own counter
% within the for loop 
m = 1;
bd = 1;
cicada_dir = cicada_home; % just for naming
num_runs = length(sub_ids);
bad_data_prefixes = ''; % initialize
samps = max([500, round(10000/num_runs)]);
for idx = 1:num_runs
    sub_id = sub_ids{idx};
    ses_id = ses_ids{idx};
    bad_data = excludes{idx}; % needs to be '1' for bad data, or '0' for OK data
    manually_adjusted = adjusteds{idx}; % needs to be '1' for manually adjsuted, '0' for not
    denoising_outlier = outliers{idx};
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
        fprintf(['   Cannot find task directory at ', task_dir, '. Skipping...\n'])
        continue;
    end

    % look for funcmask
    funcmask = [task_dir, '/funcmask.nii.gz'];
    if ~isfile(funcmask)
        fprintf(['   Cannot find funcmask at ', task_dir, '. Skipping...\n'])
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

    % a fail safe for funcmasks
    if ~isempty(dir([cleaned_dir, '/funcmask*']))
        movefile([cleaned_dir, '/funcmask*'], [cleaned_dir, '/..'])
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

    if contains(test_cleaned_file_info.name, '_CICADA_')
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

            % grab compare cleaning file now, will help find outliers
            compare_clean = readtable([task_dir, '/ic_manual_selection/compare_manual_cleaning.csv'], "ReadRowNames",true);
            % need increased GM, Smoothing Retention, best power overlap,
            % and decreased FD, DVARS, & Spikiness
            poorly_improved = (compare_clean{"GM", "Percent_Change"} < 0) | (compare_clean{"Smoothing_Retention", "Percent_Change"} < 0) | ...
                (compare_clean{"best_power_overlap_norm", "Percent_Change"} < 0) | ((compare_clean{"FD_Corr", "Percent_Change"} > 0) & ...
                (compare_clean{"DVARS_Corr", "Percent_Change"} > 0));
        else
            % else, use file_tag, it was not manually adjusted, auto cicada
            cleaned_file_info = test_cleaned_file_info;

            % grab compare cleaning file now, will help find outliers
            compare_clean = readtable([task_dir, '/ic_auto_selection/compare_auto_cleaning.csv'], "ReadRowNames",true);
            % need increased GM, Smoothing Retention, best power overlap,
            % and decreased FD, DVARS, & Spikiness
            poorly_improved = (compare_clean{"GM", "Percent_Change"} < 0) | (compare_clean{"Smoothing_Retention", "Percent_Change"} < 0) | ...
                (compare_clean{"best_power_overlap_norm", "Percent_Change"} < 0) | ((compare_clean{"FD_Corr", "Percent_Change"} > 0) & ...
                (compare_clean{"DVARS_Corr", "Percent_Change"} > 0));

        end
    else
        % this is not a cicada file, so we do not need to worry about
        % manual CICADA adjustments
        cleaned_file_info = test_cleaned_file_info;
    end

    cleaned_file = [cleaned_file_info.folder, '/', cleaned_file_info.name];
    cleaned_dir = cleaned_file_info.folder;

    % And get a helpful cleaned_file_tag 
    cleaned_file_tag = extractBetween(cleaned_file_info.name, [task_name, '_'], '_bold.nii.gz');
    cleaned_file_tag = cleaned_file_tag{:};

    % Also grab the orig file and compare file
    orig_file_info = dir([cleaned_dir, '/*_orig_*.nii.gz']);
   % make sure only one file is grabbed
    if size(orig_file_info,1) ~= 1
        fprintf('Are there no _orig_ files (or more than one) in cleaned dir?\n')
        return;
    end
    orig_file = [orig_file_info.folder, '/', orig_file_info.name];

    compare_file_info = dir([cleaned_dir, '/sub*ses*task*', task_name, '*_8p_*.nii.gz']); % need to make sure this is not the same as cleaned file
    % make sure only one file is grabbed
    if size(compare_file_info,1) ~= 1
        fprintf('Are there no _8p_ files (or more than one) in cleaned dir?\n')
        return;
    end
    compare_file = [compare_file_info.folder, '/', compare_file_info.name];

    if strcmp(compare_file, cleaned_file)
        % then the cleaned file is the 8p file... no use comparing
        % something to itself!
        compare_file = '';
    end

    if strcmp(cleaned_file, orig_file)
        % Then we are just trying to run stats on original files for
        % comparison! Do not need compare or "orig" since "cleaned" is orig
        compare_file = '';
        orig_file = '';
    end

    % Now, apply detrending and smoothing to cleaned file, orig, and compare and copy/write it to data_dir
    fprintf('Detrending & Smoothing Cleaned Data and Copying to Group Data Folder...\n')
    cleaned_file = detrend_filter_smooth(cleaned_file, funcmask, data_dir, smoothing_kernel, fpass, detrended_degree);
    if ~isempty(orig_file)
         fprintf('Detrending & Smoothing Orig Data and Copying to Group Data Folder...\n')
        orig_file = detrend_filter_smooth(orig_file, funcmask, data_dir, smoothing_kernel, fpass, detrended_degree);
    end
    % do the same to the compare_file, if it exists
    if ~isempty(compare_file)
         fprintf('Detrending & Smoothing Compare Data and Copying to Group Data Folder...\n')
        compare_file = detrend_filter_smooth(compare_file, funcmask, data_dir, smoothing_kernel, fpass, detrended_degree);
    end
    
    % OK, now FINALLY loop through cleaned_dir files, and run qc for all of them!
    % then finally individually grab relevant qc
    fprintf('Now calculating QC\n')
    [cleaned_data, data_mask, data_signal_mask, signalandnoise_overlap, qc_table, qc_corrs_table, qc_photo_paths] = cicada_get_qc(cleaned_dir, cleaned_file, samps);
    if ~isempty(orig_file)
        [~, ~, ~, ~, ~, orig_qc_corrs_table, ~] = cicada_get_qc(cleaned_dir, orig_file, samps);
    else
        orig_qc_corrs_table = table();
    end
    if ~isempty(compare_file)
        [~, ~, ~, ~, ~, compare_qc_corrs_table, ~] = cicada_get_qc(cleaned_dir, compare_file, samps);
    else
        compare_qc_corrs_table = table();
    end

    % make and get tstd
    % Also calculate and grab a tstd and grab per region
    naming_id = ['sub-', sub_id, '_ses-', ses_id, '_task-', task_name];
    call_fsl(['fslmaths ', cleaned_file, ' -Tstd -mul ', funcmask, ' ', output_dir, '/', naming_id, '_tstd.nii.gz']);
    tstd = [output_dir, '/', naming_id, '_tstd.nii.gz']; % standard deviation in time of subject. Can be helpful to look at
    
    % now, to save space, we should delete the orig and compare data from
    % data dir
    if isfile(compare_file)
        delete(compare_file)
    end
    
    if isfile(orig_file)
        delete(orig_file)
    end

    % then we need to be adding this information into the group qc stuff
    % (Group_QC struct)
    cleaned_data_list{m} = cleaned_data;
    data_mask_list{m} = data_mask;
    data_signal_mask_list{m} = data_signal_mask;
    signalandnoise_overlap_list{m} = signalandnoise_overlap;
    tstd_list{m} = tstd;
    denoising_outlier_list{m} = strcmp(denoising_outlier, '1');
    task_event_file_exist_list{m} = ~isempty(task_event_file);
    poorly_improved_list{m} = logical(poorly_improved);


    % And then move the qc photos to their respective folders & grab all qc
    % values (sampled). plot AFTER. 8p comparison is always a good base
    % comparison.
    compare_image_info = dir([task_dir, '/qc/sub*ses*task*', file_tag, '*vs*8p*qc_plots.jpg']); % specific to file tag of interest to 8p only
    compare_data_info = dir([task_dir, '/qc/sub*ses*task*', file_tag, '*vs*8p*qc_vals.mat']); % specific to file tag of interest to 8p only
    if ~isempty(compare_file)
        curr_compare_image_file = [compare_image_info.folder, '/', compare_image_info.name];
        compare_tag = extractBetween(compare_image_info.name, [task_name, '_'], '_qc_plots.jpg');
        compare_tag = compare_tag{:};
        photo_fol = [output_dir, '/qc_photos_', compare_tag];
        copyfile(curr_compare_image_file, photo_fol)
    end
    

    % Now, add all corrs into the group (compare is already done above in for loop)!
    if m == 1
         % if it is first instance, then group qc table is just subject qc
         % table
         group_qc_table = qc_table;
         group_qc_corrs_table = qc_corrs_table;
         if orig_only == 0
             group_orig_qc_corrs_table = orig_qc_corrs_table;
             group_compare_qc_corrs_table = compare_qc_corrs_table;
         end
    else
        % otherwise, tack it on to the end
        group_qc_table = [group_qc_table; qc_table];
        group_qc_corrs_table = [group_qc_corrs_table; qc_corrs_table];
        if orig_only == 0
            group_orig_qc_corrs_table = [group_orig_qc_corrs_table; orig_qc_corrs_table];
            group_compare_qc_corrs_table = [group_compare_qc_corrs_table; compare_qc_corrs_table];
        end
    end
   
    m = m+1; % increment successful counter
end
fprintf('\n Done with calculating subject data\n')


% record the bad data that was not even looked at, for easy reference later:
Group_QC.bad_data_prefixes = bad_data_prefixes;

% go to general group folder for this specific set
cd(output_dir)

% Calculate outliers based on the final qc table. Do it in three ways and
% add it to the qc table (only if this is CICADA though!):
final_qc_table = group_qc_table;
final_qc_corrs_table = group_qc_corrs_table;

% Now we can calculate outliers, but some are only calculated in this
% manner IF it is CICADA data
if cicada == 1
    % add in bold to high freq ratio, if cicada based
    boldfreq_highfreq = final_qc_table.BOLDfreq ./ final_qc_table.Highfreq;
    final_qc_table.boldfreq_highfreq_ratio = boldfreq_highfreq;
    Group_QC.boldfreq_highfreq_ratio = boldfreq_highfreq;

    % Number of total ICs (if very low -- bad data all round)
    Group_QC.low_number_total_ics = (isoutlier(final_qc_table.number_total_ics, "median")) & (final_qc_table.number_total_ics < mean(final_qc_table.number_kept_ics));
    % Low GM coverage, regular and dice 'gm_coverage_by_signal', 'signal_overlap_with_gm'
    Group_QC.low_gm_coverage_by_signal = (isoutlier(final_qc_table.gm_coverage_by_signal, "median")) & (final_qc_table.gm_coverage_by_signal < mean(final_qc_table.gm_coverage_by_signal));
    Group_QC.low_signal_overlap_with_gm = (isoutlier(final_qc_table.signal_overlap_with_gm, "median")) & (final_qc_table.signal_overlap_with_gm < mean(final_qc_table.signal_overlap_with_gm));
    Group_QC.low_gm_dice = (isoutlier(final_qc_table.gm_signal_dice, "median")) & (final_qc_table.gm_signal_dice < mean(final_qc_table.gm_signal_dice));
    
    % Fraction Signal Variance Kept (low would suggest there are very few good 
    % number of ICs (too swamped by noise) and it does not explain much of the data):
    Group_QC.low_ics_labeled_signal = (final_qc_table.number_kept_ics < 3); % If number kept ICs is less than 3, it cannot be good. And this is an extremely lenient cut off.
    Group_QC.low_fraction_signal_variance_kept = (isoutlier(final_qc_table.fraction_signal_variance_kept, "median")) & (final_qc_table.fraction_signal_variance_kept < mean(final_qc_table.fraction_signal_variance_kept));
    
    % Signal (The kept ICs should contain a fair amount of GM/Signal region
    % overlap)
    Group_QC.low_GM = (isoutlier(final_qc_table.GM, "median")) & (final_qc_table.GM < mean(final_qc_table.GM));
    
    % Smoothing
    Group_QC.low_Smoothing = (isoutlier(final_qc_table.Smoothing_Retention, "median")) & (final_qc_table.Smoothing_Retention < mean(final_qc_table.Smoothing_Retention));

    % BOLD freq / highfreq ratio
    Group_QC.low_boldfreq_highfreq_ratio = (isoutlier(boldfreq_highfreq, "median")) & (boldfreq_highfreq < mean(boldfreq_highfreq));
    
    % Power/frequency Overlap
    % if task data, take into account task power overlap
    Group_QC.low_general_power_overlap = (isoutlier(final_qc_table.general_power_overlap, "median")) & ...
        (final_qc_table.general_power_overlap < mean(final_qc_table.general_power_overlap));
    
    if sum(cell2mat(task_event_file_exist_list)) > 0 % so there are task event files!
        fprintf('Found task event files, so using best task power overlap.\n')
        Group_QC.low_besttask_power_overlap = (isoutlier(final_qc_table.besttask_power_overlap, "median")) & ...
        (final_qc_table.besttask_power_overlap < mean(final_qc_table.besttask_power_overlap));

        Group_QC.low_best_power_overlap_norm = (isoutlier(final_qc_table.best_power_overlap_norm, "median")) & ...
        (final_qc_table.best_power_overlap_norm < mean(final_qc_table.best_power_overlap_norm));

        Group_QC.low_power_overlap = Group_QC.low_general_power_overlap | Group_QC.low_boldfreq_highfreq_ratio | Group_QC.low_besttask_power_overlap;
    else
        % else it is resting state
        fprintf('No Task event files. Assumed to be resting state.\n')
        Group_QC.low_power_overlap = Group_QC.low_general_power_overlap;
    end
    
    
    % High DVARS_Corr
    Group_QC.high_DVARS_corr = (isoutlier(abs(final_qc_table.DVARS_GM_mean_abs_corr), "median")) & (abs(final_qc_table.DVARS_GM_mean_abs_corr) > median(abs(final_qc_table.DVARS_GM_mean_abs_corr)));
    
    % High FD_Corr
    Group_QC.high_FD_corr = (isoutlier(abs(final_qc_table.FD_GM_mean_abs_corr), "median")) & (abs(final_qc_table.FD_GM_mean_abs_corr) > median(abs(final_qc_table.FD_GM_mean_abs_corr)));

    % Low GM_NotGM_mean_var_prop (high GM variance and low NotGM variance
    % would occur in high signal and low noise data)
    % GM_NotGM_mean_var_prop is GM temporal variance divided by
    % NotGM temporal variance
    Group_QC.low_GM_NotGM_mean_var_prop = (isoutlier(final_qc_table.GM_NotGM_mean_var_prop, "median")) & (final_qc_table.GM_NotGM_mean_var_prop < median(final_qc_table.GM_NotGM_mean_var_prop));
    
    
    % combine to get overall outliers
    %Group_QC.cicada_outliers = logical(Group_QC.low_number_total_ics + ...
    %    Group_QC.low_fraction_signal_variance_kept + Group_QC.low_GM + ...
    %    Group_QC.low_Smoothing + Group_QC.low_power_overlap + ...
    %    Group_QC.low_gm_coverage_by_signal + Group_QC.low_signal_overlap_with_gm + Group_QC.low_gm_dice);

    % poorly improved
    Group_QC.poorly_improved = cell2mat(poorly_improved_list)';

    % Factors involving GM coverage are good, alongside factors regarding
    % power spectrum frequency (because sometimes high motion can make
    % random ICs that look like high enough GM coverage, but are not at
    % all, and are instead dominated by higher frequency noise). So, most
    % important are power overlap and GM overlap. Also anything with 2 ICs
    % as signal or less need to be cut out, as 2 is a safety number and
    % very lenient. The data must be truly horrid to not have greater than
    % 2 ICs kept as signal. Also, if the data was poorly improved (Either
    % average GM, Smoothing, or power overlap was not improved (increased), OR either
    % FD, DVARS, or spikiness was not improved (decreased).
    Group_QC.cicada_outliers = logical(Group_QC.low_gm_coverage_by_signal + Group_QC.low_gm_dice + Group_QC.low_GM_NotGM_mean_var_prop + Group_QC.low_power_overlap + Group_QC.low_boldfreq_highfreq_ratio + Group_QC.low_ics_labeled_signal + Group_QC.poorly_improved);
    

    % add to final qc table
    final_qc_table.low_number_total_ics = Group_QC.low_number_total_ics;
    final_qc_table.low_ics_labeled_signal = Group_QC.low_ics_labeled_signal; % this should give similar to low fraction signal variance kept, but this insures anything with only 2 ICs kept are excluded
    final_qc_table.low_fraction_signal_variance_kept = Group_QC.low_fraction_signal_variance_kept;
    final_qc_table.low_GM = Group_QC.low_GM; % maybe do not need if we have gm coverage by signal and signal overlap with gm
    final_qc_table.low_Smoothing = Group_QC.low_Smoothing;
    final_qc_table.low_power_overlap = Group_QC.low_power_overlap;
    final_qc_table.low_gm_dice = Group_QC.low_gm_dice;
    final_qc_table.low_gm_coverage_by_signal = Group_QC.low_gm_coverage_by_signal;
    final_qc_table.low_signal_overlap_with_gm = Group_QC.low_signal_overlap_with_gm;
    final_qc_table.low_GM_NotGM_mean_var_prop = Group_QC.low_GM_NotGM_mean_var_prop;
    final_qc_table.high_DVARS_corr = Group_QC.high_DVARS_corr;
    final_qc_table.high_FD_corr = Group_QC.high_FD_corr;
    final_qc_table.cicada_outliers = Group_QC.cicada_outliers;   
    final_qc_table.poorly_improved = Group_QC.poorly_improved;
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
writetable(final_qc_corrs_table, 'group_qc_corrs_table.csv')
Group_QC.final_qc_table = final_qc_table;
Group_QC.final_qc_corrs_table = final_qc_corrs_table;
if orig_only == 0
    % there is only a compare and orig qc corrs table if we are not just
    % looking at the original data only
    Group_QC.compare_qc_corrs_table = group_compare_qc_corrs_table;
    Group_QC.orig_qc_corrs_table = group_orig_qc_corrs_table;
else
    group_compare_qc_corrs_table = table();
    group_orig_qc_corrs_table = table();
    Group_QC.compare_qc_corrs_table = group_compare_qc_corrs_table;
    Group_QC.orig_qc_corrs_table = group_orig_qc_corrs_table;
end


% Save qc corrs data for quick and easy comparison plotting
save('qc_corrs_data.mat', 'group_qc_corrs_table', 'group_compare_qc_corrs_table', 'group_orig_qc_corrs_table' )

% Now do group plotting
dcorrt = group_qc_corrs_table; % same as final qc corrs table
ocorrt = group_orig_qc_corrs_table;

% Plotting depends on the existence of comparison tables
% Don't plot if it is just original data
% see if we have compare data to deal with here
if ~isempty(group_compare_qc_corrs_table)
    
    ccorrt = group_compare_qc_corrs_table;

    
    % Set appropriate titles and destination
    title_string = ['Group_QC_', compare_tag];
    qc_plots_dest = [output_dir, '/', title_string, '_plots.jpg'];
    
    % now we can plot
    plot_qc(cell2mat(dcorrt.Edge_Edge_Corr), cell2mat(dcorrt.FD_GM_Corr), ...
        cell2mat(dcorrt.DVARS_GM_Corr), cell2mat(dcorrt.Outbrain_Outbrain_Corr), ...
        cell2mat(dcorrt.WMCSF_WMCSF_Corr), cell2mat(dcorrt.CSF_CSF_Corr), ...
        cell2mat(dcorrt.NotGM_NotGM_Corr), cell2mat(dcorrt.GM_GM_Corr), ...
        cell2mat(dcorrt.Suscept_Suscept_Corr), cell2mat(ccorrt.Edge_Edge_Corr), ...
        cell2mat(ccorrt.FD_GM_Corr), cell2mat(ccorrt.DVARS_GM_Corr), ...
        cell2mat(ccorrt.Outbrain_Outbrain_Corr), cell2mat(ccorrt.WMCSF_WMCSF_Corr), ...
        cell2mat(ccorrt.CSF_CSF_Corr), cell2mat(ccorrt.NotGM_NotGM_Corr), ...
        cell2mat(ccorrt.GM_GM_Corr), cell2mat(ccorrt.Suscept_Suscept_Corr), ...
        cell2mat(ocorrt.Edge_Edge_Corr), cell2mat(ocorrt.FD_GM_Corr), ...
        cell2mat(ocorrt.DVARS_GM_Corr), cell2mat(ocorrt.Outbrain_Outbrain_Corr), ...
        cell2mat(ocorrt.WMCSF_WMCSF_Corr), cell2mat(ocorrt.CSF_CSF_Corr), ...
        cell2mat(ocorrt.NotGM_NotGM_Corr), cell2mat(ocorrt.GM_GM_Corr), ...
        cell2mat(ocorrt.Suscept_Suscept_Corr), ...
        [], [], [], title_string, qc_plots_dest, cleaned_file_tag)

elseif ~isempty(group_orig_qc_corrs_table)
    % we can plot now! There is no compare data
    title_string = ['Group_QC_' cleaned_file_tag, '_vs_orig'];
    qc_plots_dest = [output_dir, '/', title_string, '_plots.jpg'];
    

    plot_qc(cell2mat(dcorrt.Edge_Edge_Corr), cell2mat(dcorrt.FD_GM_Corr), ...
        cell2mat(dcorrt.DVARS_GM_Corr), cell2mat(dcorrt.Outbrain_Outbrain_Corr), ...
            cell2mat(dcorrt.WMCSF_WMCSF_Corr), cell2mat(dcorrt.CSF_CSF_Corr), ...
            cell2mat(dcorrt.NotGM_NotGM_Corr), cell2mat(dcorrt.GM_GM_Corr), ...
            cell2mat(dcorrt.Suscept_Suscept_Corr), ...
            [], [], [], [], [], [], [], [], [], ...
            cell2mat(ocorrt.Edge_Edge_Corr), cell2mat(ocorrt.FD_GM_Corr), ...
            cell2mat(ocorrt.DVARS_GM_Corr), cell2mat(ocorrt.Outbrain_Outbrain_Corr), ...
            cell2mat(ocorrt.WMCSF_WMCSF_Corr), cell2mat(ocorrt.CSF_CSF_Corr), ...
            cell2mat(ocorrt.NotGM_NotGM_Corr), cell2mat(ocorrt.GM_GM_Corr), ...
            cell2mat(ocorrt.Suscept_Suscept_Corr), ...
            [], [], [], title_string, qc_plots_dest, cleaned_file_tag)
end

fprintf('Finished Group QC Plotting!\n')


% Now work on the images:

% remove outliers from the rest of analysis
if cicada == 1
    % CICADA outliers are the outliers determined by this script. Denoising
    % outliers are the outliers that the user provided by the "outliers"
    % parameter
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

% merge temporal standard deviation files
merge_tstds_command = ['fslmerge -a ', output_dir, '/temporal_standard_deviations.nii.gz ', strjoin(tstd_list(image_keep))];
[~, ~] = call_fsl(merge_tstds_command);
% then you can delete the extras after merging!
delete([output_dir, '/sub*ses*task*tstd.nii.gz'])
[~, ~] = call_fsl(['fslmaths ', output_dir, '/temporal_standard_deviations.nii.gz -inm 1 ', output_dir '/temporal_standard_deviations.nii.gz']);

% merge func masks as well!
merge_funcmasks_command = ['fslmerge -a ', output_dir, '/funcmasks.nii.gz ', strjoin(data_mask_list(image_keep))];
[~, ~] = call_fsl(merge_funcmasks_command);

% and calculate a Tmin of data masks (where there is funcmask overlap
% across all subects)
funcmask_overlap_command = ['fslmaths ', output_dir, '/funcmasks.nii.gz -Tmin ', output_dir, '/group_funcmask.nii.gz'];
[~, ~] = call_fsl(funcmask_overlap_command);

% do the same thing, but for the data_signal_mask_list (if it exists)
if ~isempty(data_signal_mask_list)
    merge_funcmasks_signal_command = ['fslmerge -a ', output_dir, '/signal_funcmasks.nii.gz ', strjoin(data_signal_mask_list(image_keep))];
    [~, ~] = call_fsl(merge_funcmasks_signal_command);
    
    % and calculate a Tmin of data masks (where there is funcmask overlap
    % across all subects)
    funcmask_signal_overlap_command = ['fslmaths ', output_dir, '/signal_funcmasks.nii.gz -Tmin ', output_dir, '/group_signal_funcmask.nii.gz'];
    [~, ~] = call_fsl(funcmask_signal_overlap_command);
end

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
    if ~isfile('./melodic/melodic_IC.nii.gz') 
        % if melodic clearly did not finish last time, re do it.
        if isfolder('melodic')
            rmdir('melodic', 's')
        end
        fprintf(['Running: ', Group_ICA_command, '\n'])
        [~, ~] = call_fsl(Group_ICA_command);
        fprintf('Done with Group ICA (Melodic)\n\n')
    end
elseif redo_melodic == 1
    if isfolder('melodic')
        rmdir('melodic', 's')
    end

    fprintf(['Running: ', Group_ICA_command, '\n'])
    [~, ~] = call_fsl(Group_ICA_command);
    fprintf('Done with Group ICA (Melodic)\n\n')
end

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

% grab melodic file too - can be good for later images
IC_mel = niftiread('./melodic/melodic_IC.nii.gz');
IC_mel_info = niftiinfo('./melodic/melodic_IC.nii.gz');

% Update headers for 3D images
IC_3D_info = IC_probs_info;
IC_3D_info.ImageSize = IC_3D_info.ImageSize(1:end-1);
IC_3D_info.PixelDimensions = IC_3D_info.PixelDimensions(1:end-1);
IC_3D_info.Datatype = 'single';
   
% Also go through melodic output and calculate overlap % of IC to each
% network between the 7 networks, which can be useful
networks_combined = niftiread(network_file);
%networks_info = niftiinfo(network_file);

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
networks_IC_mel_max_final = networks_IC_probs_99;
for i1 = 1:size(networks,4)
    curr_dice = 0;
    curr_network_dices = dice_IC_network(:,i1);
    curr_network = networks(:,:,:,i1);
    [~, b] = sort(curr_network_dices, 'descend');

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

    % write out melodic max version, a solid reference
    curr_network_IC_mel = IC_mel(:,:,:, sort(network_ICs));
    networks_IC_mel_max_final(:,:,:,i1) = max(curr_network_IC_mel,[],4);

    networks_dice(i1) = curr_dice;
    networks_ICs{i1} = {network_ICs};
end

IC_probs_info.Datatype = 'single';
IC_probs_info.ImageSize = [IC_probs_info.ImageSize(1:end-1), size(networks,4)];
niftiwrite(cast(networks_IC_probs_99, 'single'), 'IC_probs_networks', IC_probs_info, 'Compressed', true)
niftiwrite(cast(networks_IC_probs_final, 'single'), 'IC_probs_networks_final', IC_probs_info, 'Compressed', true)
niftiwrite(cast(networks_IC_mel_max_final, 'single'), 'IC_mel_networks_final', IC_probs_info, 'Compressed', true)

Group_QC.networks_dice = networks_dice;
Group_QC.networks_ICs = networks_ICs;



save('group_qc.mat', 'Group_QC')
% When this is done, create your own table, with same image_numbers and
% image_names, with a list of what you are keeping versus not. Go through
% any you don't want to keep for qc reasons, and check if manual IC
% selection can be improved to help save it.
