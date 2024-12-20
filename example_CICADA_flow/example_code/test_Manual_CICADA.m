% to test that the manual checker works

% you will need to create an IC_manual_checker.csv that is in the same
% location, and in the same format as, IC_auto_checker.csv. The idea is
% that you open IC_auto_checker.csv, carefully review the SignalLabels (1
% is signal, 0 is noise) using manual ICA denoising guidelines, melodic
% report, and opening the melodic file on the anatomical, for example on
% fsleyes, then adjust the IC_auto_checker.csv SignalLabel(s) if needed,
% then resave the file as IC_manual_checker.csv.

% you might have to change the location of ic_manual_checker in the for
% loop depending on where you saved it, or how you named it.

CICADA_script_path = '/Users/keithdodd/GitHub/CICADA';
addpath(genpath(CICADA_script_path))
base_dir = '/Users/keithdodd/test_CICADA';


redo_mel = 0; % don't redo melodic if previously done:
compare_file = []; % will auto compare to 8p regression if given an empty array

cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
sub_ids = {'132', '135', '140'};
ses_id = '01';
task_name = 'rest';
% loop through the subjects
for i = 1:length(sub_ids)
    % grab current sub id
    sub_id = sub_ids{i};
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

    % default melodic folder location
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];

    IC_manual_checker = [output_dir, '/ic_auto_selection/IC_manual_checker.csv'];



    Manual_CICADA(output_dir, compare_file, mel_fol, IC_manual_checker)
end

