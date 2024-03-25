function [cleaned_file, data_mask, data_signal_mask, signalandnoise_overlap, qc_table, qc_corrs_table, qc_photo_paths] = cicada_get_qc(cleaned_dir, cleaned_file)
% function to grab the qc information that is needed for group qc and return relevant
% data

% may need to add to path for base scripts

if not(isfile(cleaned_file))
    fprintf(['Denoised file ', cleaned_file, 'does not exist\n'])
    return;
end

% Get cleaned dir, task dir, ses_dir, and sub_dir
cd(cleaned_dir)
cd('../')
task_dir = pwd;
[~,task_name,~]=fileparts(pwd);
cd('../') 
ses_dir = pwd;
[~,ses_id,~]=fileparts(pwd);
cd('../') 
sub_dir = pwd;
[~,sub_id,~]=fileparts(pwd);
cd(cleaned_dir)

% and get cleaned file info
cleaned_file_info = dir(cleaned_file); % mainly for the name

% figure out if cleaned file is a cicada file
if contains(cleaned_file, 'CICADA')
    cicada = 1;
else
    cicada = 0;
end

% get cicada type based on cleaned file name
if contains(cleaned_file, 'manual')
    cicada_type = 'manual';
    compare_tags = {'8p', 'auto'}; % manual will compare to 8p and auto
    adjusted = 1;
else
    cicada_type = 'auto';
    compare_tags = {'8p'}; % auto will compare only to 8p
    adjusted = 0;
end


% Grab funcmask:
if ~isfile([task_dir, '/funcmask.nii.gz'])
    fprintf(['Cannot find funcmask at ', task_dir, '/funcmask.nii.gz\n'])
    return;
end

% grab constrained funcmask if it exists too (only within brain, and not
% susceptibility low data)
if isfile([task_dir, '/funcmask_constrained.nii.gz'])
    data_mask = [task_dir, '/funcmask_constrained.nii.gz']; % use it if it is already there. This allows a user to make their own ahead of time
else
    funcfile = [task_dir, '/funcfile.nii.gz'];
    funcmask = [task_dir, '/funcmask.nii.gz'];
    anatmask = [task_dir, '/region_masks/anatmask_resam.nii.gz'];
    data_mask = make_constrained_funcmask(task_dir, funcfile, funcmask, 1);
end

% If a data_signal_mask exists, grab it (areas where signal ICs generally overlap
% with it)
data_signal_mask_path = [task_dir, '/funcmask_CICADA_', cicada_type, '_signal_constrained.nii.gz'];
if isfile(data_signal_mask_path)
    data_signal_mask = data_signal_mask_path; % constrains based on where solid signal was generally captured
else
    fprintf(['Could not find a data_signal_mask at ' data_signal_mask_path, '\n'])
    data_signal_mask = ''; % one does not exist, but this should not occur if it is CICADA
end

% grab the resampled GM mni region mask (to use for
% signalandnoise_overlap)
gm_mni_prob_info = dir([task_dir, '/region_masks/GM_prob.nii.gz']);
if size(gm_mni_prob_info,1) ~= 1
    fprintf('Cannot find gm probability file in region masks folder \n')
    return;
end
gm_mni_prob = [gm_mni_prob_info.folder, '/', gm_mni_prob_info.name];
gm_mni_thresh_array = niftiread(gm_mni_prob) > 0.67;

% Also grab helpful information from DecisionVariables*.mat
ic_select_dir = [task_dir, '/ic_', cicada_type, '_selection']; % could be auto or manual
data_vals_info = dir([ic_select_dir, '/DecisionVariables*.mat']);
data_vals = [data_vals_info.folder, '/', data_vals_info.name];
load(data_vals)
    

% Grab confounds and things that are relevant for all data:
confound_place = [task_dir, '/confounds_timeseries.csv'];
allconfounds = readtable(confound_place);
FD = table2array(allconfounds(:,{'framewise_displacement'}));
DVARS = table2array(allconfounds(:,{'dvars'}));
RMS = table2array(allconfounds(:,{'rmsd'})); 

medFD = median(FD(2:end)); % median FD to pick cut offs. Liberal is 0.55mm cut off, conservative is >0.25mm cut off
meanFD = mean(FD(2:end)); % some QC_FC likes comparing to mean FD
perc_FD_above_thresh = (sum(FD(2:end) > 0.2) / length(FD(2:end))) * 100; %percent FD greater 0.2mm, 20% cut off also used for conservative
Any_FD_above_thresh = int8(sum(FD(2:end) > 5) > 0); % also added onto conservative cut off, any FD > 5mm
medDVARS = median(DVARS(2:end));
meanRMS = mean(RMS(2:end)); % Add this into qc_group_array

% then also recalculate general correlation qc information to use here.
% Need the orig file stored in the cleaned dir, and then of course the
% actual cleaned file
orig_file_info = dir([cleaned_dir, '/*orig*.nii.gz']);
orig_file = [cleaned_dir, '/', orig_file_info.name];
[Edge_Edge_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_Outbrain_corr, WMCSF_WMCSF_corr, CSF_CSF_corr, NotGM_NotGM_corr, GM_GM_corr, Suscept_Suscept_corr, GM_mean] = CICADA_fileQC(cleaned_dir, cleaned_file, orig_file);

% grab niftiinfo to grab tr for later
cleaned_file_niftiinfo = niftiinfo(cleaned_file);
tr = cleaned_file_niftiinfo.PixelDimensions(4);

% And lets include mean and std for all relevant QC plots from last basescript:
Edge_Edge_median = median(Edge_Edge_corr);
Edge_Edge_sd = std(Edge_Edge_corr);
FD_GM_median = median(FD_GM_corr);
FD_GM_sd = std(FD_GM_corr);
DVARS_GM_median = median(DVARS_GM_corr);
DVARS_GM_sd = std(DVARS_GM_corr);
Outbrain_Outbrain_median = median(Outbrain_Outbrain_corr);
Outbrain_Outbrain_sd = std(Outbrain_Outbrain_corr);
WMCSF_WMCSF_median = median(WMCSF_WMCSF_corr);
WMCSF_WMCSF_sd = std(WMCSF_WMCSF_corr);
CSF_CSF_median = median(CSF_CSF_corr);
CSF_CSF_sd = std(CSF_CSF_corr);
NotGM_NotGM_median = median(NotGM_NotGM_corr);
NotGM_NotGM_sd = std(NotGM_NotGM_corr);
GM_GM_median = median(GM_GM_corr);
GM_GM_sd = std(GM_GM_corr);

% If it is a CICADA file, then go in and grab additional things (e.g., signal and noise ic overlap file)
if cicada == 1
    signalandnoise_overlap_info = dir([ic_select_dir, '/SignalandNoiseICOverlap.nii.gz']);
    if size(signalandnoise_overlap_info,1) ~= 1
            fprintf(['Cannot find signal and noise IC overlap file at ', [ic_select_dir, '/SignalandNoiseICOverlap.nii.gz'] '\n'])
            return;
    end
    signalandnoise_overlap = [signalandnoise_overlap_info.folder, '/', signalandnoise_overlap_info.name];

    signalandnoise_overlap_array = niftiread(signalandnoise_overlap); % will need to figure out how to adjust this given this could be TWO paths, potentially
    signal_mask_array = signalandnoise_overlap_array == 1; % positive 1 is all signal, -1 would be all noise
    gm_signal_overlap_array = gm_mni_thresh_array .* signal_mask_array;

    gm_signal_proportion_array = sum(gm_signal_overlap_array(:)) / sum(gm_mni_thresh_array(:)); % What proportion of gm is covered by signal
    signal_gm_proportion_array = sum(gm_signal_overlap_array(:)) / sum(signal_mask_array(:)); % what proportion of signal is within GM

    % and calculate a dice coefficient
    gm_and_signal_array = gm_signal_overlap_array;
    gm_or_signal_array = (gm_mni_thresh_array + signal_mask_array) > 0;
    gm_signal_dice_array = sum(gm_and_signal_array(:)) / sum(gm_or_signal_array(:));

    % grab other things only relevant for cicada files:
    percent_ICs_kept = Results.num_ICs_kept / Results.num_ICs_total; 
    percent_freq_kept_if_standardbp = 0.09 / (1 / (2* tr)); % 0.01 - 0.1 Hz
    % dof_estimate_final should be equal to percent_variance_kept *
    % percent_freq_kept * numvolumes. Percent frequency kept is
    % 100% if we did not do any frequency filtering
    DOF_estimate_if_standard_bp = Results.percent_variance_kept * percent_freq_kept_if_standardbp * Data.numvolumes;

    % Now do some variables that are only relevant for manual cicada
    if strcmp(cicada_type, 'manual')
       % also, record how many ICs were changed for manual vs auto (good to know how many were
        % "wrong") 
        changed_IC_labels_table = readtable([ic_select_dir, '/changed_signal_label_manual_ICs.csv']);
        num_ICs_adjusted = height(changed_IC_labels_table);
        percent_ICs_adjusted = (num_ICs_adjusted ./ Results.num_ICs_total) .* 100;
        percent_exp_var_adjusted = sum(Results.IC_exp_var(changed_IC_labels_table.PotentialICs)) .* 100;

    else
        % If we did manual adjustment, we should track how many ICs were
        % changed from auto. To be consistent and allow for easier coding, we
        % should also record this for auto, but it will just be 0s
        num_ICs_adjusted = 0;
        percent_ICs_adjusted = 0;
        percent_exp_var_adjusted = 0; 
    end

    % finally, qc photo paths are only relevant for cicada, and different
    % for auto vs manual
    for idx2 = 1:length(compare_tags)
         qc_photo_info = dir([task_dir, '/qc/*', cicada_type, '*cicada*', compare_tags{idx2}, '*qc_plots.jpg']);
         qc_photo_paths{idx2} = [qc_photo_info.folder, '/', qc_photo_info.name];
    end

    % Now we can combine everything in to a large table!
    % OK, now can combine all of this as a row of values

    qc_array = [{cleaned_file, cleaned_file_info.name, sub_id, ses_id, task_name, adjusted, meanRMS, ...
        medFD, meanFD, perc_FD_above_thresh, Any_FD_above_thresh, medDVARS, gm_signal_dice_array, gm_signal_proportion_array, Results.num_ICs_kept, ...
        Results.num_ICs_total, percent_ICs_kept, Results.percent_variance_kept, ...
        Data.numvolumes, DOF_estimate_if_standard_bp, ...
        NotGM_NotGM_median, NotGM_NotGM_sd, GM_GM_median, GM_GM_sd, DVARS_GM_median, DVARS_GM_sd, ...
        FD_GM_median, FD_GM_sd, CSF_CSF_median, CSF_CSF_sd, WMCSF_WMCSF_median, WMCSF_WMCSF_sd, ...
        Outbrain_Outbrain_median, Outbrain_Outbrain_sd, Edge_Edge_median, Edge_Edge_sd}, num2cell(Results.compare_cleaning.After'), ...
        {num_ICs_adjusted, percent_ICs_adjusted, percent_exp_var_adjusted}];
    
    qc_labels = ['image_path', 'image_names', 'subject', 'session', 'task', 'manually_adjusted', 'meanRMS', ...
        'median_FD', 'mean_FD', 'Percent_FD_gt_point2mm','AnyFD_gt_5mm', 'median_DVARS', 'gm_signal_dice', 'gm_coverage', 'number_kept_ics', 'number_total_ics', ...
        'fraction_kept_ics', 'fraction_signal_variance_kept', 'numvolumes', ...
        'dof_estimate_final_if_standard_resting_bp', 'NotGM_NotGM_median', 'NotGM_NotGM_sd', ...
        'GM_GM_median', 'GM_GM_sd', 'DVARS_GM_median', 'DVARS_GM_sd', 'FD_GM_median', 'FD_GM_sd', ...
        'CSF_CSF_median', 'CSF_CSF_sd', 'WMCSF_WMCSF_median', 'WMCSF_WMCSF_sd', 'Outbrain_Outbrain_median', 'Outbrain_Outbrain_sd', ...
        'Edge_Edge_median', 'Edge_Edge_sd', Results.compare_cleaning.Properties.RowNames', 'Num_ICs_Adjusted', 'Percent_ICs_Adjusted', 'Percent_Exp_Var_Adjusted'];
    
    qc_table = cell2table(qc_array, 'VariableNames', qc_labels);

else
    % if not cicada, the following variables should be blank:
    signalandnoise_overlap = '';
    qc_photo_paths = '';

    % Now we can combine everything in to a large table!
    % OK, now can combine all of this as a row of values

    qc_array = {cleaned_file, cleaned_file_info.name, sub_id, ses_id, task_name, meanRMS, ...
        medFD, meanFD, perc_FD_above_thresh, Any_FD_above_thresh, medDVARS, ...
        Data.numvolumes, NotGM_NotGM_median, NotGM_NotGM_sd, GM_GM_median, GM_GM_sd, DVARS_GM_median, DVARS_GM_sd, ...
        FD_GM_median, FD_GM_sd, CSF_CSF_median, CSF_CSF_sd, WMCSF_WMCSF_median, WMCSF_WMCSF_sd, ...
        Outbrain_Outbrain_median, Outbrain_Outbrain_sd, Edge_Edge_median, Edge_Edge_sd};
    
    qc_labels = ["image_path", "image_names", "subject", "session", "task", "meanRMS", ...
        "median_FD", "mean_FD", "Percent_FD_gt_point2mm","AnyFD_gt_5mm", "median_DVARS", ...
        "numvolumes", "NotGM_NotGM_median", "NotGM_NotGM_sd", "GM_GM_median", "GM_GM_sd", ...
        "DVARS_GM_median", "DVARS_GM_sd", "FD_GM_median", "FD_GM_sd", ...
        "CSF_CSF_median", "CSF_CSF_sd", "WMCSF_WMCSF_median", "WMCSF_WMCSF_sd", ...
        "Outbrain_Outbrain_median", "Outbrain_Outbrain_sd", "Edge_Edge_median", "Edge_Edge_sd"];
    
    qc_table = cell2table(qc_array, 'VariableNames', qc_labels);
end


% Grab a sampling of the correlation values, because otherwise it might
% become too much
samps = 500;
qc_corrs_table = grab_corr_sampling(Edge_Edge_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_Outbrain_corr, WMCSF_WMCSF_corr, CSF_CSF_corr, NotGM_NotGM_corr, GM_GM_corr, Suscept_Suscept_corr, samps);

end