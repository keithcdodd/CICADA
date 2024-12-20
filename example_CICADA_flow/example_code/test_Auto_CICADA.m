% try running Automatic CICADA (calls Auto_CICADA) on the test data to see if things are
% seemingly configured and working properly!

% Note: After Auto_CICADA is successful you have two choices:
% (1) Complete Automatic CICADA testing by running test_Group_CICADA.m
% which runs cicada_group_qc.m to finish with Group CICADA. Make sure
% file_tag = '_auto_'
% (2) If you want to fully perform manual ICA denoising through CICADA, you can do the following: 
% (2a) Manually adjust the SignalLabels in IC_auto_checker.csv in the ic_auto_checker directory for each cicada dataset and save the changes as IC_manual_checker.csv,
% (2b) Then run test_Manual_CICADA.m
% (2c) Now just run Group CICADA with test_Group_CICADA.m. Make sure
% file_tag = '_manual_'; 

% Note, after this runs successfully, you can also test Group CICADA
% through the test_Group_CICADA.m which calls cicada_group_qc.m

% Calling either fmriprep_auto_CICADA.m or Auto_CICADA.m should work with
% the appropriate inputs if everything is configured correctly

% Auto_CICADA.m is more applicable to different groupings of data, so will
% use that here

% test Auto_CICADA.m on the resting state example dataset!
% Note: You will have to make sure the base directory is the test_CICADA
% directory, and that the path to the CICADA scripts is correct (and/or
% added)

CICADA_script_path = '/Users/keithdodd/GitHub/CICADA';
addpath(genpath(CICADA_script_path))
base_dir = '/Users/keithdodd/test_CICADA';


redo_mel = 0; % don't redo melodic if previously done:
mel_fol = [];


cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
sub_ids = {'132', '135', '140'};
ses_id = '01';
task_name = 'rest';
fmriprep_dir = [base_dir, '/bids_data/derivatives/fmriprep'];

% loop through the subjects
for i = 1:length(sub_ids)
    % grab current sub id
    sub_id = sub_ids{i};
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

    % default melodic folder location
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    
    % grab functionals
    func_dir = [fmriprep_dir, '/sub-', sub_id, '/ses-', ses_id, '/func'];
    funcfile = [func_dir, '/sub-', sub_id, '_ses-01_task-rest_space-MNI152NLin6Asym_res-02_desc-preproc_bold.nii.gz'];
    funcmask = [func_dir, '/sub-', sub_id, '_ses-01_task-rest_space-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz'];
    confoundsfile = [func_dir, '/sub-', sub_id, '_ses-01_task-rest_desc-confounds_timeseries.tsv'];

    compare_file = []; % will auto compare to 8p regression if given an empty array
    task_events_file = []; % does not apply to resting state
    
    % grab anatomicals
    anat_dir = [fmriprep_dir, '/sub-', sub_id, '/ses-', ses_id, '/anat'];
    anatfile = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_desc-preproc_T1w.nii.gz'];
    anatmask = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz'];
    gm_prob = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_label-GM_probseg.nii.gz'];
    wm_prob = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_label-WM_probseg.nii.gz'];
    csf_prob = [anat_dir, '/sub-', sub_id, '_ses-01_space-MNI152NLin6Asym_res-02_label-CSF_probseg.nii.gz'];

    % default tolerance
    tolerance = 5;

    Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, redo_mel, mel_fol, compare_file, task_events_file, anatfile, anatmask, gm_prob, wm_prob, csf_prob, tolerance)
end
