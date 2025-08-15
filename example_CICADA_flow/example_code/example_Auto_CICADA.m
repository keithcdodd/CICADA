% Example code that would work with the files and file structure shown
% in the example image. This is provided to help user(s) set up their own
% data and function calls.

% Note: After Auto_CICADA is successful you have two choices:
% (1) Complete Automatic CICADA testing by running something similar to example_Group_CICADA.m
% which runs cicada_group_qc.m to finish with Group CICADA. Make sure
% file_tag = '_auto_'
% (2) If you want to fully perform manual ICA denoising through CICADA, you can do the following: 
% (2a) Manually adjust the SignalLabels in IC_auto_checker.csv in the ic_auto_checker directory for each cicada dataset and save the changes as IC_manual_checker.csv,
% (2b) Then run example_Manual_CICADA.m
% (2c) Now just run Group CICADA with example_Group_CICADA.m. Make sure
% file_tag = '_manual_'; 

% Calling either fmriprep_auto_CICADA.m or Auto_CICADA.m should work with
% the appropriate inputs if everything is configured correctly

% Auto_CICADA.m is more applicable to different groupings of data, so will
% use that here


CICADA_script_path = '/Users/keithdodd/GitHub/CICADA';
addpath(genpath(CICADA_script_path))
base_dir = '/Users/keithdodd/test_CICADA';


redo_mel = 0; % don't redo melodic if previously done:
mel_fol = [];


cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
sub_ids = {'132', '135', '140'};

% to compare to your own data, you can set compare_files to a different
% full filepath for each subject you want to compare (same length as
% sub_ids)
% alternatively, you can feed it one of the valid compare_tags for standard
% comparison: '6p', '8p', '12p', '16p', '18p', '24p', '28p', '30p', '32p', or '36p'
% '8p' is the standard 8 parameter regression (6 motion + CSF + WM)!
compare_files = {'8p'};

% For task events file, can either supply one fullfile path that applies to
% all data being examined, or one unique task event file per subject/data
% (similar to compare_files)
task_events_files = {};

ses_id = '01';
task_name = 'rest';
fmriprep_dir = [base_dir, '/bids_data/derivatives/fmriprep'];

% Set despike to 1 if you want to lightly despike the fMRI data before running ICA
% Note: the idea is to very lightly denoise the fMRI data to potentially
% assist with IC decomposition. This is largely untested, but theoretically
% may help with particularly noisy data and should have low threshold to
% hurt the analysis. This was not used as part of the original paper.
despike = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(compare_files)
    compare_files = {'8p'}; % if just doing default or forgot to set the parameter
end

if isempty(task_events_files)
    task_events_files = []; % no task event files being used
end

if length(compare_files) ~= length(sub_ids)
    if length(compare_files) ~= 1
        % assuming you are running one session and one task at a time here
        fprintf('Compare Files must either be the same length as the subjects you are looking at, or length of 1\n')
        return;
    end
end

if length(task_events_files) ~= length(sub_ids)
    if length(task_events_files) ~= 1 || ~isempty(task_events_files)
        % assuming you are running one session and one task at a time here
        fprintf('Task event files must either be the same length as the subjects you are looking at, or length of 1, or length of 0\n')
        return;
    end
end

% loop through the subjects
for i = 1:length(sub_ids)
    % grab current sub id
    sub_id = sub_ids{i};
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

    % comparison(s)
    if length(compare_files) > 1
        compare_file = compare_files{i};
    elseif isempty(compare_files)
        compare_file = '8p';
    else
        compare_file = compare_files{1};
    end

    % task event(s)
    if length(task_events_files) > 1
        task_events_file = task_events_files{i};
    elseif isempty(task_events_files)
        task_events_file = [];
    else
        task_events_file = task_events_files{1};
    end

    % default melodic folder location
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    
    % grab functionals
    func_dir = [fmriprep_dir, '/sub-', sub_id, '/ses-', ses_id, '/func'];
    funcfile = [func_dir, '/sub-', sub_id, '_ses-01_task-rest_space-MNI152NLin6Asym_res-02_desc-preproc_bold.nii.gz'];
    funcmask = [func_dir, '/sub-', sub_id, '_ses-01_task-rest_space-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz'];
    confoundsfile = [func_dir, '/sub-', sub_id, '_ses-01_task-rest_desc-confounds_timeseries.tsv'];
    
    % grab anatomicals
    anat_dir = [fmriprep_dir, '/sub-', sub_id, '/ses-', ses_id, '/anat'];
    anatfile = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_desc-preproc_T1w.nii.gz'];
    anatmask = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz'];
    gm_prob = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_label-GM_probseg.nii.gz'];
    wm_prob = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_label-WM_probseg.nii.gz'];
    csf_prob = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_label-CSF_probseg.nii.gz'];

	if despike == 1
        % lightly despike the functional data first!
        fprintf('Lightly Despiking the Data Before IC Decomposition!\n')
        robust_z_thresh = 4;
        [funcfile_despiked] = despike_fMRI(funcfile, gm_prob, robust_z_thresh);
        funcfile = funcfile_despiked;
    end
    
    % default tolerance
    tolerance = 5;

    Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, redo_mel, mel_fol, compare_file, task_events_file, anatfile, anatmask, gm_prob, wm_prob, csf_prob, tolerance)
end
