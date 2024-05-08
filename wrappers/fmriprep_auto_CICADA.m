function fmriprep_auto_CICADA(fmriprep_dir, cicada_dir, sub_id, ses_id, task_name, anat_ses_id, redo_mel, mel_fol, task_events_file, compare_file, tolerance)
% A wrapper script to make it easier to work with fmriprep datasets and run
% CICADA on them. 
% fmriprep_dir: home directory for fmriprep folder: e.g. /path/fmriprep
% cicada_dir: home directory for cicada folder e.g., /path/cicada
% sub_id: e.g., '101'
% ses_id e.g., '01'
% task_name e.g., 'visual_run-01' (include run tags in this if they exist!)
% anat_ses_id: Whichever ses_id  has the best anat folder to use
% e.g., '01'
% compare file: What you want the QC plots to compare CICADA to... default
% is 8param denoising
% task_events_file: Which task events file to use for this, should be a
% .tsv, in the bids format specified in bids/fmriprep (basically columns
% for onset, duration, and trial_type is most common).
% mel_fol: MELODIC folder to use, usually with CICADA it will run MELODIC
% for you, but if you already have a MELODIC folder you insist on using,
% you can supply that here


% check if fmriprep has sessions or not (if only one session, might not!)
if isfolder([fmriprep_dir, '/sub-', sub_id, '/ses-', ses_id])
    has_ses = 1; % this is the expected structure, with multiple sessions
else
    has_ses = 0; 
    ses_id = ''; % then give empty ses_id
    anat_ses_id = '';
end

% Now check for non-necessary variables
if ~exist('redo_mel', 'var') || (redo_mel ~= 0 && redo_mel ~=1)
    redo_mel = 0; % default is to not redo melodic
end

if ~exist('compare_file', 'var') || ~ischar(compare_file) || isempty(compare_file)
    compare_file=[];
    compare_file_record = 'Standard 8 parameter compare';
elseif ~isfile(compare_file)
    fprintf(['Compare file not found at ', compare_file, '\n'])
    return;
else
    compare_file_record = compare_file;
end

% if task_events_file does not exist, make it empty array [], 'x' would
% work too given how I wrote the scripts.
if ~exist('task_events_file', 'var') || ~ischar(task_events_file) || isempty(task_events_file)
    task_events_file=[]; 
    task_events_file_record = 'None Provided';
elseif ~isfile(task_events_file)
    fprintf(['Task events file not found at ', task_events_file, '\n'])
    return;
else
    task_events_file_record = task_events_file;
end

if ~exist('mel_fol', 'var') || ~ischar(mel_fol) || isempty(mel_fol)
    if has_ses == 1
        mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic']; % will look in default location
    else
        % if no designated session, just make it session 01
        mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-01/', task_name, '/melodic']; % will look in default location
    end
elseif ~isfolder(mel_fol)
    fprintf(['MELODIC folder not found at ', mel_fol, '\n'])
    if has_ses == 1
        mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    else
        % if no designated session, just make it session 01
        mel_fol = [cicada_dir, '/sub-', sub_id, '/ses-01/', task_name, '/melodic'];
    end
    fprintf(['Will run new melodic in new default folder: ' mel_fol, '\n'])
end



% First do all the checks to make sure everything seems legit!
% check that expected fmriprep directories exist first
if has_ses == 1
    fmriprep_func_dir = [fmriprep_dir, '/sub-', sub_id, '/ses-', ses_id, '/func'];
    fmriprep_anat_dir = [fmriprep_dir, '/sub-', sub_id, '/ses-', anat_ses_id, '/anat'];
else
    % If there is no session folder in fmriprep, work with that
    fmriprep_func_dir = [fmriprep_dir, '/sub-', sub_id, '/func'];
    fmriprep_anat_dir = [fmriprep_dir, '/sub-', sub_id, '/anat'];
end


if ~isfolder(fmriprep_func_dir)
    fprintf(['Cannot find fmriprep func dir at ', fmriprep_func_dir, '\n'])
    return;
elseif ~isfolder(fmriprep_anat_dir)
    fprintf(['Cannot find fmriprep anat dir at ', fmriprep_anat_dir, '\n'])
    return;
else
    % first grab things from func folder
    cd(fmriprep_func_dir)
    funcfile_info = dir(['*', sub_id, '*-', ses_id, '*', task_name, '*', 'space-MNI*preproc_bold.nii.gz']);
    funcfile = [funcfile_info.folder, '/', funcfile_info.name];
    if ~isfile(funcfile)
        fprintf(['Cannot find fmriprep funcfile at ', funcfile, '\n'])
        return;
    end

    funcmask_info = dir(['*', sub_id, '*-', ses_id, '*', task_name, '*', 'space-MNI*brain_mask.nii.gz']);
    funcmask = [funcmask_info.folder, '/', funcmask_info.name];
    confounds_info =  dir(['*', sub_id, '*-', ses_id, '*', task_name, '*', 'desc-confounds_timeseries.tsv']);
    confoundsfile = [confounds_info.folder, '/', confounds_info.name];

    if ~isfile(confoundsfile)
        fprintf(['Cannot find confounds_file at ', confoundsfile, '\n'])
        return;
    end
    % and now grab from best anat folder
    % Note, if you use a different structural other than T1w, you may need
    % to change that below
    cd(fmriprep_anat_dir)
    anatfile_info = dir(['*', sub_id, '*-', anat_ses_id, '*', 'space-MNI*preproc_T1w.nii.gz']);
    anatfile = [anatfile_info.folder, '/', anatfile_info.name];
    if ~isfile(funcfile)
        fprintf(['Cannot find fmriprep anatfile at ', anatfile, '\n'])
        return;
    end
    anatmask_info = dir(['*', sub_id, '*-', anat_ses_id, '*', 'space-MNI*desc-brain_mask.nii.gz']);
    anatmask = [anatmask_info.folder, '/', anatmask_info.name];
    gm_prob_info = dir(['*', sub_id, '*-', anat_ses_id, '*', 'space-MNI*label-GM_probseg.nii.gz']);
    gm_prob = [gm_prob_info.folder, '/', gm_prob_info.name];
    wm_prob_info = dir(['*', sub_id, '*-', anat_ses_id, '*', 'space-MNI*label-WM_probseg.nii.gz']);
    wm_prob = [wm_prob_info.folder, '/', wm_prob_info.name];
    csf_prob_info = dir(['*', sub_id, '*-', anat_ses_id, '*', 'space-MNI*label-CSF_probseg.nii.gz']);
    csf_prob = [csf_prob_info.folder, '/', csf_prob_info.name];
end

if has_ses == 1
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];
else
    % just put in a ses-01 folder for consistency
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-01/', task_name];
end

if ~exist('tolerance', 'var')
    tolerance = 3;
end

% Output everything to check:
fprintf(['\nfmriprep_dir: ', fmriprep_dir, '\n'])
fprintf(['cicada_dir: ', cicada_dir, '\n'])
fprintf(['sub_id: ', sub_id, '\n'])
fprintf(['ses_id: ', ses_id, '\n'])
fprintf(['task_name: ', task_name, '\n'])
fprintf(['ses_id for anat: ', anat_ses_id, '\n'])
fprintf(['task_events_file: ', task_events_file_record, '\n'])
fprintf(['compare_file: ', compare_file_record, '\n'])
fprintf(['melodic folder: ', mel_fol, '\n\n'])


% For melodic folder, we can just check if the melodic folder exists, does
% it have one of the final inputs?

Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, redo_mel, mel_fol, compare_file, task_events_file, anatfile, anatmask, gm_prob, wm_prob, csf_prob, tolerance)

end