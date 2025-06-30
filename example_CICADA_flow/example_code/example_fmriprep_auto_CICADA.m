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
sub_ids = {'107', '158'};
ses_id = '01'; % which scan has the functional scan
anat_ses_id = '01'; % which session has the anatomical scan of interest
task_name = 'rest';
fmriprep_dir = [base_dir, '/bids_data/derivatives/fmriprep'];

% to compare to your own data, you can set compare_files to a different
% full filepath for each subject you want to compare (same length as
% sub_ids)
% alternatively, you can feed it one of the valid compare_tags for standard
% comparison: '6p', '8p', '12p', '16p', '18p', '24p', '32p', or '36p'
% '8p' is the standard 8 parameter regression (6 motion + CSF + WM)!
compare_files = {'8p'};

% For task events file, can either supply one fullfile path that applies to
% all data being examined, or one unique task event file per subject/data
% (similar to compare_files)
task_events_files = {};

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
    if length(task_events_files) ~= 1 && ~isempty(task_events_files)
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

    % task event(s)
    if length(task_events_files) > 1
        task_events_file = task_events_files{i};
    elseif isempty(task_events_files)
        task_events_file = [];
    else
        task_events_file = task_events_files{1};
    end
    
    % default tolerance
    tolerance = 5;

    fmriprep_auto_CICADA(fmriprep_dir, cicada_dir, sub_id, ses_id, task_name, anat_ses_id, redo_mel, mel_fol, task_events_file, compare_file, tolerance)
end

