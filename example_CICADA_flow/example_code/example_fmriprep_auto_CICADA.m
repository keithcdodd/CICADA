% Example code that would work with the files and file structure shown
% in the example image. This is provided to help user(s) set up their own
% data and function calls. This function uses the file structure to
% properly call Auto_CICADA wrapper function. 


CICADA_script_path = '/Users/keithdodd/GitHub/CICADA';
addpath(genpath(CICADA_script_path))
base_dir = '/Users/keithdodd/test_CICADA';


redo_mel = 0; % don't redo melodic if previously done:
mel_fol = [];


cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
sub_ids = {'132', '135', '140'};
ses_id = '01'; % which scan has the functional scan
anat_ses_id = '01'; % which session has the anatomical scan of interest
task_name = 'rest';
fmriprep_dir = [base_dir, '/bids_data/derivatives/fmriprep'];

% loop through the subjects
for i = 1:length(sub_ids)
    % grab current sub id
    sub_id = sub_ids{i};
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

    % default melodic folder location
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    
    compare_file = []; % will auto compare to 8p regression if given an empty array
    task_events_file = []; % does not apply to resting state
    
    % default tolerance
    tolerance = 5;

    fmriprep_auto_CICADA(fmriprep_dir, cicada_dir, sub_id, ses_id, task_name, anat_ses_id, redo_mel, mel_fol, task_events_file, compare_file, tolerance)
end

