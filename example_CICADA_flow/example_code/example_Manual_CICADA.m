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

cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
sub_ids = {'132', '135', '140'};
ses_id = '01';
task_name = 'rest';

% to compare to your own data, you can set compare_files to a different
% full filepath for each subject you want to compare (same length as
% sub_ids)
% alternatively, you can feed it one of the valid compare_tags for standard
% comparison: '6p', '8p', '12p', '16p', '18p', '24p', '32p', or '36p'
% '8p' is the standard 8 parameter regression (6 motion + CSF + WM)!
compare_files = {'8p'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(compare_files)
    compare_files = {'8p'}; % if just doing default or forgot to set the parameter
end

if length(compare_files) ~= length(sub_ids)
    if length(compare_files) ~= 1
        % assuming you are running one session and one task at a time here
        fprintf('Compare Files must either be the same length as the subjects you are looking at, or length of 1\n')
        return;
    end
end


% loop through the subjects
for i = 1:length(sub_ids)
    % grab current sub id
    sub_id = sub_ids{i};
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

    % default melodic folder location
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];

    % comparison(s)
    if length(compare_files) > 1
        compare_file = compare_files{i};
    elseif isempty(compare_files)
        compare_file = '8p'; % defaults to 8p comparison
    else
        compare_file = compare_files{1};
    end

    IC_manual_checker = [output_dir, '/ic_auto_selection/IC_manual_checker.csv'];

    Manual_CICADA(output_dir, compare_file, mel_fol, IC_manual_checker)
end

