function [cleaned_data, data_mask, signalandnoise_overlap, qc_table, qc_corrs_table, qc_photo_paths] = cicada_get_qc(cleaned_file)
% function to grab the qc information that is needed for group qc and return relevant
% data

% may need to add to path for base scripts

if not(isfile(cleaned_file))
    fprintf(['Denoised file ', cleaned_file, 'does not exist\n'])
    return;
end

cleaned_file_info = dir(cleaned_file);
cleaned_dir = cleaned_file_info.folder;

% Get cleaned dir, task dir, ses_dir, and subj_dir
cd(cleaned_dir)
cd('../')
task_dir = pwd;
cd('../') 
ses_dir = pwd;
cd('../') 
subj_dir = pwd;
cd(cleaned_dir)

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

data_mask = [task_dir, '/funcmask.nii.gz'];

% grab the resampled GM mni region mask (to use for
% signalandnoise_overlap)
gm_mni_prob_info = dir([task_dir, '/region_masks/GM_prob.nii.gz']);
if size(gm_mni_prob_info,1) ~= 1
    fprintf('Cannot find gm probability file inr region masks folder \n')
    return;
end
gm_mni_prob = [gm_mni_prob_info.folder, '/', gm_mni_prob_info.name];
gm_mni_thresh_array = niftiread(gm_mni_prob) > 0.67;

% Also grab helpful information from DecisionVariables*.mat
ic_select_dir = [task_dir, '/ic_', cicada_type, '_selection']; % could be auto or manual
data_vals = dir([ic_select_dir, '/DecisionVariables*.mat']);
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
orig_file_info = dir([cleaned_dir, '*orig*.nii.gz']);
orig_file = [cleaned_dir, '/', orig_file_info.name];
[Edge_GM_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, NotGM_GM_corr, GM_GM_autocorr, GM_mean] = CICADA_fileQC(cleaned_file, orig_file);

% And lets include mean and std for all relevant QC plots from last basescript:
Edge_GM_median = median(Edge_GM_corr);
Edge_GM_sd = std(Edge_GM_corr);
FD_GM_median = median(FD_GM_corr);
FD_GM_sd = std(FD_GM_corr);
DVARS_GM_median = median(DVARS_GM_corr);
DVARS_GM_sd = std(DVARS_GM_corr);
Outbrain_GM_median = median(Outbrain_GM_corr);
Outbrain_GM_sd = std(Outbrain_GM_corr);
WMCSF_GM_median = median(WMCSF_GM_corr);
WMCSF_GM_sd = std(WMCSF_GM_corr);
CSF_GM_median = median(CSF_GM_corr);
CSF_GM_sd = std(CSF_GM_corr);
NotGM_GM_median = median(NotGM_GM_corr);
NotGM_GM_sd = std(NotGM_GM_corr);
GM_GM_median = median(GM_GM_autocorr);
GM_GM_sd = std(GM_GM_autocorr);

% If it is a CICADA file, then go in and grab additional things (e.g., signal and noise ic overlap file)
if cicada == 1
    signalandnoise_overlap_info = dir([task_dir, '/', ic_select_dir, '/SignalandNoiseICOverlap.nii.gz']);
    if size(signalandnoise_overlap_info,1) ~= 1
            fprintf('Cannot find signal and noise IC overlap file \n')
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

    qc_array = [{cleaned_file_info.name, subj_id, ses_id, task_name, adjusted, meanRMS, ...
        medFD, meanFD, perc_FD_above_thresh, Any_FD_above_thresh, medDVARS, gm_signal_dice_array, gm_signal_proportion_array, Results.num_ICs_kept, ...
        Results.num_ICs_total, percent_ICs_kept, Results.percent_variance_kept, ...
        Data.numvolumes, DOF_estimate_final, DOF_estimate_if_standard_bp, ...
        NotGM_GM_median, NotGM_GM_sd, GM_GM_median, GM_GM_sd, DVARS_GM_median, DVARS_GM_sd, ...
        FD_GM_median, FD_GM_sd, CSF_GM_median, CSF_GM_sd, WMCSF_GM_median, WMCSF_GM_sd, ...
        Outbrain_GM_median, Outbrain_GM_sd, Edge_GM_median, Edge_GM_sd}, num2cell(Results.compare_cleaning.After'), ...
        {num_ICs_adjusted, percent_ICs_adjusted, percent_exp_var_adjusted}];
    
    qc_labels = ['image_names', 'subject', 'session', 'task', 'manually_adjusted', 'meanRMS', ...
        'median_FD', 'mean_FD', 'Percent_FD_gt_point2mm','AnyFD_gt_5mm', 'median_DVARS', 'gm_signal_dice', 'gm_coverage', 'number_kept_ics', 'number_total_ics', ...
        'fraction_kept_ics', 'fraction_signal_variance_kept', 'numvolumes', ...
        'dof_estimate_final', 'dof_estimate_final_if_standard_resting_bp', 'NotGM_GM_median', 'NotGM_GM_sd', ...
        'GM_GM_median', 'GM_GM_sd', 'DVARS_GMmedian', 'DVARS_GMsd', 'FD_GMmedian', 'FD_GMsd', ...
        'CSF_GMmedian', 'CSF_GMsd', 'WMCSF_GMmedian', 'WMCSF_GMsd', 'Outbrain_GMmedian', 'Outbrain_GMsd', ...
        'Edge_GMmedian', 'Edge_GMsd', Results.compare_cleaning.Properties.RowNames', 'Num_ICs_Adjusted', 'Percent_ICs_Adjusted', 'Percent_Exp_Var_Adjusted'];
    
    qc_table = cell2table(qc_array, 'VariableNames', qc_labels);

else
    % if not cicada, the following variables should be blank:
    signalandnoise_overlap = '';
    qc_photo_paths = '';

    % Now we can combine everything in to a large table!
    % OK, now can combine all of this as a row of values

    qc_array = {cleaned_file_info.name, subj_id, ses_id, task_name, meanRMS, ...
        medFD, meanFD, perc_FD_above_thresh, Any_FD_above_thresh, medDVARS, ...
        Data.numvolumes, NotGM_GM_median, NotGM_GM_sd, GM_GM_median, GM_GM_sd, DVARS_GM_median, DVARS_GM_sd, ...
        FD_GM_median, FD_GM_sd, CSF_GM_median, CSF_GM_sd, WMCSF_GM_median, WMCSF_GM_sd, ...
        Outbrain_GM_median, Outbrain_GM_sd, Edge_GM_median, Edge_GM_sd};
    
    qc_labels = ['image_names', 'subject', 'session', 'task', 'meanRMS', ...
        'median_FD', 'mean_FD', 'Percent_FD_gt_point2mm','AnyFD_gt_5mm', 'median_DVARS', ...
        'numvolumes', 'NotGM_GM_median', 'NotGM_GM_sd', 'GM_GM_median', 'GM_GM_sd', ...
        'DVARS_GM_median', 'DVARS_GM_sd', 'FD_GM_median', 'FD_GM_sd', ...
        'CSF_GM_median', 'CSF_GM_sd', 'WMCSF_GM_median', 'WMCSF_GM_sd', ...
        'Outbrain_GM_median', 'Outbrain_GM_sd', 'Edge_GM_median', 'Edge_GM_sd'];
    
    qc_table = cell2table(qc_array, 'VariableNames', qc_labels);
end


% Grab a sampling of the correlation values, because otherwise it might
% become too much
samps=500;
NotGM_GM_randperm = randperm(size(NotGM_GM_corr, 1)); % and then grab ~samps of these
NotGM_GM_randperm = NotGM_GM_randperm(1:samps);

GM_GM_randperm = randperm(size(GM_GM_autocorr, 1)); % and then grab ~samps of these
GM_GM_randperm = GM_GM_randperm(1:samps);

DVARS_GMrandperm = randperm(size(DVARS_GM_corr, 1)); % and then grab ~samps of these
DVARS_GMrandperm = DVARS_GMrandperm(1:samps);

FD_GMrandperm = randperm(size(FD_GM_corr, 1)); % and then grab ~samps of these
FD_GMrandperm = FD_GMrandperm(1:samps);

CSF_GMrandperm = randperm(size(CSF_GM_corr, 1)); % and then grab ~samps of these
CSF_GMrandperm = CSF_GMrandperm(1:samps);

WMCSF_GMrandperm = randperm(size(WMCSF_GM_corr, 1)); % and then grab ~samps of these
WMCSF_GMrandperm = WMCSF_GMrandperm(1:samps);

Outbrain_GMrandperm = randperm(size(Outbrain_GM_corr, 1)); % and then grab ~samps of these
Outbrain_GMrandperm = Outbrain_GMrandperm(1:samps);

Edge_GMrandperm = randperm(size(Edge_GM_corr, 1)); % and then grab ~samps of these
Edge_GMrandperm = Edge_GMrandperm(1:samps);

% We will make a qc_corrs_table that holds manual (if it exists or is
% calculated, otherwise 0s), auto, 8p, & 9p, because why not?

% Now that we have the randomized samps numbers, we can actually grab them
% from the correlations
NotGM_GM_corr = NotGM_GM_corr(NotGM_GM_randperm);
GM_GM_corr = GM_GM_corr(GM_GM_randperm);
DVARS_GM_corr = DVARS_GM_corr(DVARS_GMrandperm);
FD_GM_corr = FD_GM_corr(FD_GMrandperm);
CSF_GM_corr = CSF_GM_corr(CSF_GMrandperm);
WMCSF_GM_corr = WMCSF_GM_corr(WMCSF_GMrandperm);
Outbrain_GM_corr = Outbrain_GM_corr(Outbrain_GMrandperm);
Edge_GM_corr = Edge_GM_corr(Edge_GMrandperm);

qc_corrs_array = [Edge_GM_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, NotGM_GM_corr, GM_GM_corr];
qc_corrs_labels = {'Edge_GM_Corr', 'FD_GM_Corr', ...
    'DVARS_GM_Corr', 'Outbrain_GM_Corr', 'WMCSF_GM_Corr', ...
    'CSF_GM_Corr', 'NotGM_GM_Corr', 'GM_GM_AutoCorr'};
qc_corrs_table = array2table(qc_corrs_array, 'VariableNames',qc_corrs_labels);

end