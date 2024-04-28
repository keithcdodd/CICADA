function fmriprep_manual_CICADA(cicada_dir, sub_id, ses_id, task_name, IC_manual_checker, compare_file, mel_fol)
% A wrapper script to make it easier to work with fmriprep datasets and run
% CICADA on them. 
% cicada_dir: home directory for cicada folder e.g., /path/cicada
% sub_id: e.g., '101'
% ses_id e.g., '01'
% task_name e.g., 'visual_run-01' (include run tags in this if they exist!)
% IC_manual_checker: filepath to your updated IC_manual_checker.csv where
% you manipulated and renamed the IC_auto_checker.csv to do manual IC
% denoising.
% compare file: What you want the QC plots to compare CICADA to... default
% is 8param denoising
% mel_fol: MELODIC folder to use, usually with CICADA it will run MELODIC
% for you, but if you already have a MELODIC folder you insist on using,
% you can supply that here

if isempty(ses_id)
    % then the original data did not have sessions, so default to session
    % 01 to match other outputs
    ses_id = '01';
end

% Now check for non-necessary variables
if ~exist('compare_file', 'var') || ~ischar(compare_file) || isempty(compare_file)
    compare_file=[];
    compare_file_record = 'Standard 8 parameter compare';
elseif ~isfile(compare_file)
    fprintf(['Compare file not found at ', compare_file, '\n'])
    return;
else
    compare_file_record = compare_file;
end

if ~exist('mel_fol', 'var') || ~ischar(mel_fol) || isempty(mel_fol)
    mel_fol=[]; % will end up using default
elseif ~isfolder(mel_fol)
    fprintf(['MELODIC folder not found at ', mel_fol, '\n'])
    mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    fprintf(['Will run new melodic in new default folder: ' mel_fol, '\n'])
end

output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];

if ~isfolder(output_dir)
    fprintf(['Cannot find output_dir (task_dir) at ' output_dir, '. Is your folder structure wrong perhaps?'])
    return;
end

if ~exist('IC_manual_checker', 'var') || isempty(IC_manual_checker) || ~isfile(IC_manual_checker)
    % If not given assume it is in its normal spot and look for it
    if isfile([output_dir, '/ic_auto_selection/IC_manual_checker.csv'])
        IC_manual_checker = [output_dir, '/ic_auto_selection/IC_manual_checker.csv'];
    elseif isfile([output_dir, '/ic_manual_selection/IC_manual_checker.csv'])
            IC_manual_checker = [output_dir, '/ic_manual_selection/IC_manual_checker.csv'];
    else
        fprintf('Cannot find IC_manual_checker. Exiting...\n')
        return
    end
end

% Output everything to check:
fprintf(['cicada_dir: ', cicada_dir, '\n'])
fprintf(['sub_id: ', sub_id, '\n'])
fprintf(['ses_id: ', ses_id, '\n'])
fprintf(['task_name: ', task_name, '\n'])
fprintf(['compare_file: ', compare_file_record, '\n'])
fprintf(['melodic folder: ', mel_fol, '\n\n'])


% For melodic folder, we can just check if the melodic folder exists, does
% it have one of the final inputs?

Manual_CICADA(output_dir, compare_file, mel_fol, IC_manual_checker)

end