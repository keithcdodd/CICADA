function cicada_group_qc(cicada_home, group_qc_home, task_name, output_dirname, file_tag, voxelwise_scale_flag, smoothing_kernel, fpass, detrended_degree, redo_melodic, sub_ids, ses_ids, excludes, adjusteds, compare_tag, task_event_files)
% function to run group qc
% cicada_home: is the cicada home directory, the general home input folder
% group_qc_home: the general group qc output dir
% output_dirname: what name do you want for your output directory? This
% will end up like the following: group_qc_home/task_name/output_dirname
% task_name: The task id for the group qc data, e.g. 'visual_run-01', or
% 'rest'
% sub_ids, ses_ids, excludes, adjusteds, and
% task_event_files should all be cell arrays of characters of the same number of
% rows. This function will loop through each row of them all and match them
% for the call. 
% file_tag: needs to be unique to the type of denoised/cleaned file you are
% putting together in group qc. e.g., '_auto_', '_manual_', '_8p_'
% compare_rag: similar to file tag, for the file you are comparing against,
% e.g., '8p' being the default

% voxelwise_scale_flag should be a 1 or a 0. If 1, voxelwise scaling is 
% performed (after detrending, bandpassing, and smoothing). 
% Otherwise, voxelwise scaling is not performed 

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
% adjusteds (either 0 or 1 for manually IC adjusted) e.g., {'0', '1', '0'}
% task_event_files should contain paths to the task_event_file(s)
% (optional), also as a cell array

% if no task_event_file_list is given (e.g., it is all resting state) then
% create an empty cell array of chars the same length as sub_id_list
if ~exist('task_event_files', 'var') || isempty(task_event_files)
    task_event_files = cell(size(sub_ids));
    task_event_files(:) = {''}; 
end

if ~exist('compare_tag', 'var') || isempty(compare_tag)
    compare_tag = '8p'; % just use standard 8p comparison if it is not provided by user
end

if isfile(compare_tag)
    fprintf('Error. You seem to have provided a whole singular compare file, but it needs to be an identifiable tag of a file in the cleaned dir!\n')
    return
end

valid_tags = {'6p', '8p', '9p', '12p', '16p', '18p', '24p', '28p', '30p', '32p', '36p'};
if ~ismember(compare_tag, valid_tags)
    compare_tag = '8p'; % set to default if it is not an appropriate tag
end

% now compare_tag is an appropriate one!
fprintf(['Comparing to ', compare_tag, '.\n'])

% now check that all lists are the same length, and that all are lists of
% char arrays
if ~isequal(length(sub_ids), length(ses_ids), length(excludes), length(adjusteds), length(task_event_files))
    fprintf('Your lists (sub_ids, ses_ids, excludes, adjusteds, task_event_files) exist at different lengths... quitting\n')
    return
end

test_class = {''};

if ~isequal(class(test_class), class(sub_ids), class(ses_ids), class(excludes), class(adjusteds), class(task_event_files))
    fprintf('Your lists are not all cell arrays \n')
    return
end

if ~isequal(class(test_class{1}), class(sub_ids{1}), class(ses_ids{1}), class(excludes{1}), class(adjusteds{1}), class(task_event_files{1}))
    fprintf('Your lists do not contain all chars (should be cell arrays of characters!) \n')
    return
end

% if no voxelwise_scale_flag exists, make it 0 (no voxelwise_scaling)
if ~exist('voxelwise_scale_flag', 'var') || ~isnumeric(voxelwise_scale_flag)
    voxelwise_scale_flag = 0;
end

% default to no voxelwise scaling, only a 1 will result in doing it. Better
% safe than sorry!
if voxelwise_scale_flag ~= 1
    fprintf('Default no voxelwise scaling will be applied.\n')
    voxelwise_scale_flag = 0;
end

% if no smoothing kernel exists, make it default value
default_smoothing = 0;
if ~exist('smoothing_kernel', 'var') || ~isnumeric(smoothing_kernel) || smoothing_kernel < 0
    fprintf('Default smoothing gaussian kernel (FWHMx = 1.5x voxel size) will be applied\n')
    default_smoothing = 1; % flag to make it the default later on in the function
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
    fpass = []; % if bandpass frequencies are not readable, assume default of no bandpassing
end

% Make separate photo folders, so it is easier to just scroll through them
% all in the future. There is only 8p compare if auto, but 8p and auto if
% manual 

% initialize keeping track of data issues
data_notes_columnNames = {'Subject', 'Session', 'Task', 'Type', 'Note'};
data_notes_types = {'string', 'string', 'string', 'string', 'string'};
data_notes = table('Size', [0, length(data_notes_columnNames)], 'VariableNames', data_notes_columnNames, 'VariableTypes', data_notes_types);
data_notes_iter = 1;

% go into first participant's info to extract the correct naming, just compare to 8p here for ease:
first_task_dir = [cicada_home, '/sub-', sub_ids{1}, '/ses-', ses_ids{1}, '/', task_name];
compare_image_info = dir([first_task_dir, '/qc/sub*ses*task*', file_tag, '*vs*', compare_tag, '*qc_plots.jpg']); % specific to file tag of interest, can be multiples


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
    comparison_tag = extractBetween(compare_image_info(h).name, [task_name, '_'], '_qc_plots.jpg');
    comparison_tag = comparison_tag{:};
    photo_fol = [output_dir, '/qc_photos_', comparison_tag];
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
voxel_size = round(mean(first_image_nifi_info.PixelDimensions(1:3)));

if default_smoothing == 1
   smoothing_kernel = 1.5 * voxel_size; % default smoothing is 1.5 times the voxel size
   fprintf(['Will smooth at FWHM of 1.5x voxel size: ', num2str(smoothing_kernel) ,' mm.\n'])
end

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
    task_event_file = task_event_files{idx};

    fprintf(['\nRunning sub ', sub_id, ' ses ', ses_id, ' task ', task_name, '\n'])

    % if bad data, whole run needs to be skipped and not included in qc
    % (use continue)
    if ~strcmp(bad_data, '0')
        % this is bad data from the get go (recognized previously, like
        % with mriqc, or at time of scan), do not include in group qc
        % analysis
        bad_data_prefixes{bd} = ['sub-', sub_id, '_ses-', ses_id, '_task-', task_name];
        data_notes{data_notes_iter, :} = {sub_id, ses_id, task_name, 'Excluded', 'Marked to Exclude by User Before Denoising'};
        data_notes_iter = data_notes_iter + 1; % increment for table
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
            data_notes{data_notes_iter, :} = {sub_id, ses_id, task_name, 'Manually Adjusted', 'Manual IC classification performed by User, with manual denoising used instead of automated because of better denoising performance.'};
            data_notes_iter = data_notes_iter + 1; % increment for table

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

    compare_file_info = dir([cleaned_dir, '/sub*ses*task*', task_name, '*_' compare_tag, '_*.nii.gz']); % need to make sure this is not the same as cleaned file
    % make sure only one file is grabbed
    if size(compare_file_info,1) ~= 1
        fprintf(['Are there no _', compare_tag, '_ files (or more than one) in cleaned dir? Check naming conventions!\n'])
        return;
    end
    compare_file = [compare_file_info.folder, '/', compare_file_info.name];

    % if strcmp(compare_file, cleaned_file)
    %     % then the cleaned file is the 8p file... no use comparing
    %     % something to itself!
    %     compare_file = '';
    % end
    % 
    % if strcmp(cleaned_file, orig_file)
    %     % Then we are just trying to run stats on original files for
    %     % comparison! Do not need compare or "orig" since "cleaned" is orig
    %     compare_file = '';
    %     orig_file = '';
    % end

    % Now, apply detrending and smoothing to cleaned file, orig, and compare and copy/write it to data_dir
    % smoothing kernel is FWHM mm, as this is then properly converted
    % within the detrend_filter_smooth function
    fprintf('Detrending & Smoothing Cleaned Data and Copying to Group Data Folder...\n')
    cleaned_file = detrend_filter_smooth(cleaned_file, funcmask, data_dir, smoothing_kernel, fpass, detrended_degree);
    
    fprintf('Detrending & Smoothing Compare Data and Copying to Group Data Folder...\n')
    compare_file = detrend_filter_smooth(compare_file, funcmask, data_dir, smoothing_kernel, fpass, detrended_degree);

    fprintf('Detrending & Smoothing Orig Data and Copying to Group Data Folder...\n')
    orig_file = detrend_filter_smooth(orig_file, funcmask, data_dir, smoothing_kernel, fpass, detrended_degree);

    if voxelwise_scale_flag == 1
        fprintf('Performing Voxel-wise Scaling to Mean 100\n')
        cleaned_file = voxelwise_scale(cleaned_file); % voxelwise scale each voxel to mean 100
        compare_file = voxelwise_scale(compare_file); % voxelwise scale each voxel to mean 100
        orig_file = voxelwise_scale(orig_file); % voxelwise scale each voxel to mean 100
    end
    
    
    % OK, now FINALLY loop through cleaned_dir files, and run qc for all of them!
    % then finally individually grab relevant qc
    fprintf('Now calculating QC\n')
    [cleaned_data, data_mask, data_signal_mask, signalandnoise_overlap, cleaned_qc_table, ...
        cleaned_qc_corrs_table, cleaned_qc_photo_paths, ...
        compare_qc_table, compare_qc_corrs_table, ...
        orig_qc_table, orig_qc_corrs_table] = cicada_get_qc(cleaned_dir, cleaned_file, compare_file, orig_file, samps);
    
    % if ~isempty(orig_file)
    %     [~, ~, ~, ~, ~, orig_qc_corrs_table, ~] = cicada_get_qc(cleaned_dir, orig_file, samps);
    % else
    %     orig_qc_corrs_table = table();
    % end
    % if ~isempty(compare_file)
    %     [~, ~, ~, ~, ~, compare_qc_corrs_table, ~] = cicada_get_qc(cleaned_dir, compare_file, samps);
    % else
    %     compare_qc_corrs_table = table();
    % end

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
    task_event_file_exist_list{m} = ~isempty(task_event_file);
    poorly_improved_list{m} = logical(poorly_improved);


    % QC info!
    % And then move the qc photos to their respective folders & grab all qc
    % values (sampled). plot AFTER. 8p comparison is always a good base
    % comparison.
    compare_image_info = dir([task_dir, '/qc/sub*ses*task*', file_tag, '*vs*' compare_tag, '*qc_plots.jpg']); % specific to file tag of interest to 8p only
    compare_data_info = dir([task_dir, '/qc/sub*ses*task*', file_tag, '*vs*', compare_tag, '*qc_vals.mat']); % specific to file tag of interest to 8p only
    if ~isempty(compare_file)
        curr_compare_image_file = [compare_image_info.folder, '/', compare_image_info.name];
        comparison_tag = extractBetween(compare_image_info.name, [task_name, '_'], '_qc_plots.jpg');
        comparison_tag = comparison_tag{:};
        photo_fol = [output_dir, '/qc_photos_', comparison_tag];
        copyfile(curr_compare_image_file, photo_fol)
    end

    % now is a good time to grab network identifiability too!
    network_identifiability_file = [task_dir, '/qc/network_identifiability.nii.gz'];
    if isfile(network_identifiability_file)
        fprintf('Copying network identifiability image information:\n')
        curr_network_identifiability = niftiread(network_identifiability_file);
        curr_denoised_network_identifiability = curr_network_identifiability(:,:,:,1);
        curr_compare_network_identifiability = curr_network_identifiability(:,:,:,2);
        curr_orig_network_identifiability = curr_network_identifiability(:,:,:,3);

        % save organized data. Will need to write to nifti once complete
        denoised_network_identifiability(:,:,:,m) = curr_denoised_network_identifiability;
        compare_network_identifiability(:,:,:,m) = curr_compare_network_identifiability;
        orig_network_identifiability(:,:,:,m) = curr_orig_network_identifiability;
    else
        fprintf('No network identifiability image found!\n')
    end

    % Now, add all corrs into the group (compare is already done above in for loop)!
    if m == 1
         % if it is first instance, then group qc table is just subject qc
         % table
         group_qc_table = cleaned_qc_table;
         group_qc_corrs_table = cleaned_qc_corrs_table;
         if orig_only == 0
             group_compare_qc_table = compare_qc_table;
             group_orig_qc_table = orig_qc_table;

             group_orig_qc_corrs_table = orig_qc_corrs_table;
             group_compare_qc_corrs_table = compare_qc_corrs_table;
         end
    else
        % otherwise, tack it on to the end
        group_qc_table = [group_qc_table; cleaned_qc_table];
        group_qc_corrs_table = [group_qc_corrs_table; cleaned_qc_corrs_table];
        if orig_only == 0
            group_orig_qc_table = [group_orig_qc_table; orig_qc_table];
            group_compare_qc_table = [group_compare_qc_table; compare_qc_table];

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

    % if the mean is above 0.15 in both, it often looks bad in QC!
    Group_QC.high_retained_motion = (abs(final_qc_table.DVARS_GM_mean_abs_corr) >= 0.15) | (abs(final_qc_table.FD_GM_mean_abs_corr) >= 0.15); % 0.15 is when QC plots tend to really break down visually. This is still liberally keeping data.


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
    cicada_outliers = logical(Group_QC.high_retained_motion + Group_QC.low_gm_dice + Group_QC.low_GM_NotGM_mean_var_prop + Group_QC.low_power_overlap + Group_QC.low_boldfreq_highfreq_ratio + Group_QC.low_ics_labeled_signal + Group_QC.poorly_improved);
    Group_QC.cicada_outliers = cicada_outliers;

    % record into data notes:
    for ico = 1:length(cicada_outliers)
        if cicada_outliers(ico)
            data_notes{data_notes_iter, :} = {sub_id, ses_id, task_name, 'CICADA Outlier', 'Group CICADA automatically labeled this as an outlier compared to the other data.'};
            data_notes_iter = data_notes_iter + 1; % increment for table
        end
    end

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
    final_qc_table.high_retained_motion = Group_QC.high_retained_motion;
    final_qc_table.cicada_outliers = Group_QC.cicada_outliers;   
    final_qc_table.poorly_improved = Group_QC.poorly_improved;
end

% Now, get liberal outliers (not restrictive)
Group_QC.liberal_outliers = final_qc_table.mean_FD > 0.55;

% Finally, conservative outliers: median FD
Group_QC.conservative_outliers = (final_qc_table.mean_FD > 0.25) | ...
    (final_qc_table.Percent_FD_gt_point2mm > 20) | ...
    (final_qc_table.AnyFD_gt_5mm);

% sort data_notes
data_notes_sorted = sortrows(data_notes, 'Subject', 'ascend');
% save data notes, and write it as a table!
Group_QC.data_notes = data_notes_sorted;
writetable(data_notes_sorted, 'group_qc_data_notes_table.csv')

% Make a final list of what images (out of the ones kept in group_qc_data)
% are used for future analyses
% remove outliers from the rest of analysis
if cicada == 1
    % CICADA outliers are the outliers determined by this script. 
    image_keep = ~Group_QC.cicada_outliers;
else
    image_keep = true(size(Group_QC.liberal_outliers)); % keep everyone if you do not have cicada outliers
end

% save image_keep so it is easy in the future to know what images were left
% out of MELODIC Group IC
Group_QC.Images_Used = image_keep;

% add to final qc table
final_qc_table.liberal_outliers = Group_QC.liberal_outliers;
final_qc_table.conservative_outliers = Group_QC.conservative_outliers;
final_qc_table.image_keep = Group_QC.Images_Used; % All images except cicada_outliers

writetable(final_qc_table, 'group_qc_table.csv') % can use this and sort by gm_signal and signal_gm proportions and also number_kept_ics to help find the bad or adjustable data
writetable(final_qc_corrs_table, 'group_qc_corrs_table.csv') % does have subject, session, and task recorded as column variables too!
Group_QC.final_qc_table = final_qc_table;
Group_QC.final_qc_corrs_table = final_qc_corrs_table;
if orig_only == 0
    % there is only a compare and orig qc corrs table if we are not just
    % looking at the original data only
    Group_QC.compare_qc_table = group_compare_qc_table;
    Group_QC.orig_qc_table = group_orig_qc_table;

    Group_QC.compare_qc_corrs_table = group_compare_qc_corrs_table;
    Group_QC.orig_qc_corrs_table = group_orig_qc_corrs_table;
else
    group_compare_qc_table = table();
    group_orig_qc_table = table();

    group_compare_qc_corrs_table = table();
    group_orig_qc_corrs_table = table();

    Group_QC.compare_qc_table = group_compare_qc_table;
    Group_QC.orig_qc_table = group_orig_qc_table;

    Group_QC.compare_qc_corrs_table = group_compare_qc_corrs_table;
    Group_QC.orig_qc_corrs_table = group_orig_qc_corrs_table;
end

% also make a qc_corrs table that is without the cicada outliers!

% Step 1: Identify bad subject-session-task combinations
cicada_outlier_combos = group_qc_table(image_keep == 0, {'subject', 'session', 'task'}); % Marked by CICADA as a "cicada_outlier" 

% Step 2: Start with full table, then remove bad rows
group_qc_corrs_cicada_outliers_removed_table = group_qc_corrs_table;
group_compare_qc_corrs_cicada_outliers_removed_table = group_compare_qc_corrs_table;
group_orig_qc_corrs_cicada_outliers_removed_table = group_orig_qc_corrs_table;

% Step 3: Loop through each bad combination and remove matches
for i = 1:height(cicada_outlier_combos)
    is_bad_cicada_cleaned = strcmp(group_qc_corrs_cicada_outliers_removed_table.subject, cicada_outlier_combos.subject{i}) & ...
             strcmp(group_qc_corrs_cicada_outliers_removed_table.session, cicada_outlier_combos.session{i}) & ...
             strcmp(group_qc_corrs_cicada_outliers_removed_table.task, cicada_outlier_combos.task{i});

    is_bad_cicada_compare = strcmp(group_compare_qc_corrs_cicada_outliers_removed_table.subject, cicada_outlier_combos.subject{i}) & ...
             strcmp(group_compare_qc_corrs_cicada_outliers_removed_table.session, cicada_outlier_combos.session{i}) & ...
             strcmp(group_compare_qc_corrs_cicada_outliers_removed_table.task, cicada_outlier_combos.task{i});

    is_bad_cicada_orig = strcmp(group_orig_qc_corrs_cicada_outliers_removed_table.subject, cicada_outlier_combos.subject{i}) & ...
             strcmp(group_orig_qc_corrs_cicada_outliers_removed_table.session, cicada_outlier_combos.session{i}) & ...
             strcmp(group_orig_qc_corrs_cicada_outliers_removed_table.task, cicada_outlier_combos.task{i});

    group_qc_corrs_cicada_outliers_removed_table(is_bad_cicada_cleaned, :) = [];  % remove matching rows
    group_compare_qc_corrs_cicada_outliers_removed_table(is_bad_cicada_compare, :) = [];  % remove matching rows
    group_orig_qc_corrs_cicada_outliers_removed_table(is_bad_cicada_orig, :) = [];  % remove matching rows
    
end


% Save qc corrs data for quick and easy comparison plotting
save('qc_corrs_data.mat', 'group_qc_corrs_table', 'group_compare_qc_corrs_table', 'group_orig_qc_corrs_table', ...
    'group_qc_corrs_cicada_outliers_removed_table', 'group_compare_qc_corrs_cicada_outliers_removed_table', 'group_orig_qc_corrs_cicada_outliers_removed_table')
save('qc_data.mat', 'group_qc_table', 'group_compare_qc_table', 'group_orig_qc_table')

% Now do group plotting, first with all data
dcorrt = group_qc_corrs_table; % same as final qc corrs table
ocorrt = group_orig_qc_corrs_table;

dcorrt_cor = group_qc_corrs_cicada_outliers_removed_table; % cicada outliers removed
ocorrt_cor = group_orig_qc_corrs_cicada_outliers_removed_table;


% Plotting depends on the existence of comparison tables
% Don't plot if it is just original data
% see if we have compare data to deal with here
if ~isempty(group_compare_qc_corrs_table)
    
    ccorrt = group_compare_qc_corrs_table;
    ccorrt_cor = group_compare_qc_corrs_cicada_outliers_removed_table; % cicada_outliers removed

    
    % Set appropriate titles and destination
    title_string = ['Group_QC_', comparison_tag];
    qc_plots_dest = [output_dir, '/', title_string, '_plots.jpg'];
    
    % now we can plot everyone
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

   
    % Set appropriate titles and destination
    title_string_cor = ['Group_QC_', comparison_tag, '_cicada_outliers_removed'];
    qc_plots_dest_cor = [output_dir, '/', title_string_cor, '_plots.jpg'];
    
    % now we can plot with cicada_outliers removed for comparison
    plot_qc(cell2mat(dcorrt_cor.Edge_Edge_Corr), cell2mat(dcorrt_cor.FD_GM_Corr), ...
        cell2mat(dcorrt_cor.DVARS_GM_Corr), cell2mat(dcorrt_cor.Outbrain_Outbrain_Corr), ...
        cell2mat(dcorrt_cor.WMCSF_WMCSF_Corr), cell2mat(dcorrt_cor.CSF_CSF_Corr), ...
        cell2mat(dcorrt_cor.NotGM_NotGM_Corr), cell2mat(dcorrt_cor.GM_GM_Corr), ...
        cell2mat(dcorrt_cor.Suscept_Suscept_Corr), cell2mat(ccorrt_cor.Edge_Edge_Corr), ...
        cell2mat(ccorrt_cor.FD_GM_Corr), cell2mat(ccorrt_cor.DVARS_GM_Corr), ...
        cell2mat(ccorrt_cor.Outbrain_Outbrain_Corr), cell2mat(ccorrt_cor.WMCSF_WMCSF_Corr), ...
        cell2mat(ccorrt_cor.CSF_CSF_Corr), cell2mat(ccorrt_cor.NotGM_NotGM_Corr), ...
        cell2mat(ccorrt_cor.GM_GM_Corr), cell2mat(ccorrt_cor.Suscept_Suscept_Corr), ...
        cell2mat(ocorrt_cor.Edge_Edge_Corr), cell2mat(ocorrt_cor.FD_GM_Corr), ...
        cell2mat(ocorrt_cor.DVARS_GM_Corr), cell2mat(ocorrt_cor.Outbrain_Outbrain_Corr), ...
        cell2mat(ocorrt_cor.WMCSF_WMCSF_Corr), cell2mat(ocorrt_cor.CSF_CSF_Corr), ...
        cell2mat(ocorrt_cor.NotGM_NotGM_Corr), cell2mat(ocorrt_cor.GM_GM_Corr), ...
        cell2mat(ocorrt_cor.Suscept_Suscept_Corr), ...
        [], [], [], title_string_cor, qc_plots_dest_cor, cleaned_file_tag)

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

    % Set appropriate titles and destination
    title_string_cor = ['Group_QC_', cleaned_file_tag, '_vs_orig_cicada_outliers_removed'];
    qc_plots_dest_cor = [output_dir, '/', title_string_cor, '_plots.jpg'];
    
    % now we can plot with cicada_outliers removed for comparison
    plot_qc(cell2mat(dcorrt_cor.Edge_Edge_Corr), cell2mat(dcorrt_cor.FD_GM_Corr), ...
        cell2mat(dcorrt_cor.DVARS_GM_Corr), cell2mat(dcorrt_cor.Outbrain_Outbrain_Corr), ...
        cell2mat(dcorrt_cor.WMCSF_WMCSF_Corr), cell2mat(dcorrt_cor.CSF_CSF_Corr), ...
        cell2mat(dcorrt_cor.NotGM_NotGM_Corr), cell2mat(dcorrt_cor.GM_GM_Corr), ...
        cell2mat(dcorrt_cor.Suscept_Suscept_Corr), ...
        [], [], [], [], [], [], [], [], [], ...
        cell2mat(ocorrt_cor.Edge_Edge_Corr), cell2mat(ocorrt_cor.FD_GM_Corr), ...
        cell2mat(ocorrt_cor.DVARS_GM_Corr), cell2mat(ocorrt_cor.Outbrain_Outbrain_Corr), ...
        cell2mat(ocorrt_cor.WMCSF_WMCSF_Corr), cell2mat(ocorrt_cor.CSF_CSF_Corr), ...
        cell2mat(ocorrt_cor.NotGM_NotGM_Corr), cell2mat(ocorrt_cor.GM_GM_Corr), ...
        cell2mat(ocorrt_cor.Suscept_Suscept_Corr), ...
        [], [], [], title_string_cor, qc_plots_dest_cor, cleaned_file_tag)
end

fprintf('Finished Group QC Plotting!\n')

% Now write network identifiability nifti files!
if isfile(network_identifiability_file)
    network_identifiability_info = niftiinfo(network_identifiability_file);
    network_identifiability_info.ImageSize = [network_identifiability_info.ImageSize(1:3), (m-1)];
    network_identifiability_info.Datatype = 'single';
    
    % denoised
    network_identifiability_info.Filename = [output_dir, '/network_identifiability_', cleaned_file_tag, '.nii.gz'];
    niftiwrite(single(denoised_network_identifiability), network_identifiability_info.Filename, network_identifiability_info, 'Compressed', true);
    
    % compare
    compare_file_tag = extractAfter(compare_tag, '_vs_');
    network_identifiability_info.Filename = [output_dir, '/network_identifiability_', compare_file_tag, '.nii.gz'];
    niftiwrite(single(compare_network_identifiability), network_identifiability_info.Filename, network_identifiability_info, 'Compressed', true);
    
    % orig
    network_identifiability_info.Filename = [output_dir, '/network_identifiability_orig.nii.gz'];
    niftiwrite(single(orig_network_identifiability), network_identifiability_info.Filename, network_identifiability_info, 'Compressed', true);
end


% Now work on the images (we already have image_keep set from earlier)

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

% Potential Future Idea: After group MELODIC, make automatic decisions on signal or
% noise like ICs, based on spatial map and back projections of power
% spectrum.. etc.

% Run Group MELODIC
% First, make text file of the image names
% will need a text file where separated by commas (no spaces) and then in
% the melodic call $(cat )
% Make sure to not include data files that were designated as bad enough to
% remove
data_files = cleaned_data_list;
all_data_files = data_files;
data_files = all_data_files(image_keep); % removing cicada_outliers if they exist
writecell(data_files, 'image_names.txt')
group_funcmask_file = fullfile(pwd, 'group_funcmask.nii.gz');

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
IC_probs_mask = logical(IC_probs > 0.99); % keeps only very best overlap

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

Group_QC.networks = networks;
Group_QC.IC_probs_mask = IC_probs_mask; % 99% mask
Group_QC.data_files = data_files;
network_names = {'Visual', 'Sensorimotor', 'Dorsal Attention', 'Salience', 'Executive', 'Default', 'Subcortical'};
IC_assignment_table = assign_ICs_to_networks_clust(networks, IC_probs_mask, network_names);

% Also just write out the melodic ICs that are maintained overall!
kept_ICs_mel = IC_mel(:,:,:, IC_assignment_table.IC_Index);
IC_mel_info.Datatype = 'single';
IC_mel_info.ImageSize = [IC_mel_info.ImageSize(1:end-1), length(IC_assignment_table.IC_Index)];
niftiwrite(cast(kept_ICs_mel, 'single'), 'IC_mel_networks', IC_mel_info, 'Compressed', true)

% record table
writetable(IC_assignment_table, 'IC_assignment_table.csv', 'Delimiter',',')

networks_IC_probs = IC_probs(:,:,:,IC_assignment_table.IC_Index);
niftiwrite(cast(networks_IC_probs, 'single'), 'IC_probs_networks', IC_mel_info, 'Compressed', true) % max of probability files for related networks


% make files for masking each network and maximum z score too
networks_IC_probs_mask = zeros(size(networks));
networks_IC_mel_max = zeros(size(networks));
for i1 = 1:max(IC_assignment_table.Assigned_Network)

    curr_network_ICs = IC_assignment_table.IC_Index(IC_assignment_table.Assigned_Network == i1);

    % writeout IC prob network file:
    networks_IC_probs_mask(:,:,:,i1) = max(IC_probs(:,:,:, curr_network_ICs),[], 4);

    % write out melodic max version, a solid reference
    curr_network_IC_mel = IC_mel(:,:,:, curr_network_ICs);
    [~, max_pos_abs_z] = max(abs(curr_network_IC_mel),[],4, 'linear');
    networks_IC_mel_max(:,:,:,i1) = curr_network_IC_mel(max_pos_abs_z); % will follow the sign now
end

IC_3D_info = IC_probs_info;
IC_3D_info.Datatype = 'single';
IC_3D_info.ImageSize = [IC_probs_info.ImageSize(1:end-1), size(networks,4)];

niftiwrite(cast(networks_IC_probs_mask, 'single'), 'IC_mask_networks', IC_3D_info, 'Compressed', true) % mask of max 99%+ probability for related networks
niftiwrite(cast(networks_IC_mel_max, 'single'), 'IC_mel_networks_zmax', IC_3D_info, 'Compressed', true) % maximum (positive or negative) z values for ICs in network


% OK, and then we can also just sort group IC from most to least gray
% matter
% Grab DVARS and FD though so you can easily measure correlations too
subj_func_files = data_files;
DVARS_all = {};
FD_all = {};
notGM_all = {};
for i1 = 1:length(data_files)
    curr_subj_home = [cicada_home, '/', group_qc_table.subject{i1}, '/', group_qc_table.session{i1}, '/', group_qc_table.task{i1}];

    % grab notGM mask
    curr_notGM_mask = [curr_subj_home, '/region_masks/NotGM_mask.nii.gz'];
    notGM_all{end+1} = curr_notGM_mask;

    if exist([curr_subj_home, '/confounds_timeseries.tsv'], 'file') ~= 0
        confound_place = [curr_subj_home, '/confounds_timeseries.tsv'];
        % if it is a .tsv, then it should be tab delimited.
        allconfounds = readtable(confound_place, 'FileType', 'text', 'Delimiter', '\t');
    elseif exist([curr_subj_home, '/confounds_timeseries.csv'], 'file') ~= 0
        confound_place = [curr_subj_home, '/confounds_timeseries.csv'];
        % matlab can natively read csv no problem. CSV is probably better.
        allconfounds = readtable(confound_place);
    else
        fprintf('ERROR: Cannot find confounds_timeseries as either csv or tsv?\n')
        return;
    end
    curr_DVARS = table2array(allconfounds(:,{'dvars'}));
    curr_FD = table2array(allconfounds(:,{'framewise_displacement'}));

    DVARS_all{end+1} = curr_DVARS;
    FD_all{end+1} = curr_FD;

    %% grab HRF information:
    mat_file = [curr_subj_home, '/ic_auto_selection/DecisionVariables_Auto.mat'];
    load(mat_file)

    % Only load HRF_task if it exists
    
    if isfield(Data, 'HRF_task')
        hrf_general_power = Data.HRF_general.hrf_power_interp_norm;
        hrf_task_power = Data.HRF_task.hrf_power_norm;
        hrf_conditions = Data.Task.hrf_conditions;
    else
        % only general power!

        hrf_general_power = Data.HRF_general.hrf_power_interp_norm;
        hrf_task_power = []; % just set it to be the same thing for ease (keep in mind, will normally have multiple columns)
        hrf_conditions = [];
    end
end

% Now output general group IC results with information that can be helpful
% in deciding if true signal or not. (Power spectra, GM overlap, timeseries
% correlations)
% make sure to test with task data too
T_results = rank_ICs_by_group_NSP(output_dir, IC_mel, group_funcmask_file, subj_func_files, notGM_all, tr, FD_all, DVARS_all, hrf_general_power, hrf_task_power, hrf_conditions);

save('group_qc.mat', 'Group_QC', 'T_results')

