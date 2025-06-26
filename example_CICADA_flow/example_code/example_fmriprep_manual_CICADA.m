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

    fmriprep_manual_CICADA(cicada_dir, sub_id, ses_id, task_name, IC_manual_checker, compare_file, mel_fol)
end

