% Example code that would work with the files and file structure shown
% in the example image. This is provided to help user(s) set up their own
% data and function calls. This function uses the file structure to
% properly call Manual_CICADA wrapper function. 

% may need to change location/name for IC_manual_checker in for loop


CICADA_script_path = '/Users/keithdodd/GitHub/CICADA';
addpath(genpath(CICADA_script_path))
base_dir = '/Users/keithdodd/test_CICADA';

mel_fol = [];

cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
sub_ids = {'132', '135', '140'};
ses_id = '01'; % which scan has the functional scan
task_name = 'rest';

% loop through the subjects
for i = 1:length(sub_ids)
    % grab current sub id
    sub_id = sub_ids{i};
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

    % default melodic folder location
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    
    compare_file = []; % will auto compare to 8p regression if given an empty array

    IC_manual_checker = [output_dir, '/ic_auto_selection/IC_manual_checker.csv'];

    fmriprep_manual_CICADA(cicada_dir, sub_id, ses_id, task_name, IC_manual_checker, compare_file, mel_fol)
end

