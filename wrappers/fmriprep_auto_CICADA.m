function fmriprep_auto_CICADA(fmriprep_dir, cicada_dir, subj_id, ses_id, task_name, anat_ses_id, task_events_file, compare_file, mel_fol)
% A wrapper script to make it easier to work with fmriprep datasets and run
% CICADA on them. 
% fmriprep_dir: home directory for fmriprep folder: e.g. /path/fmriprep
% cicada_dir: home directory for cicada folder e.g., /path/cicada
% subj_id: e.g., '101'
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

% Now check for non-necessary variables
if ~exist('compare_file', 'var') || ~ischar(compare_file) || isempty(compare_file)
    compare_file=[];
elseif ~isfile(compare_file)
    fprintf(['Compare file not found at ', compare_file, '\n'])
    return
end

% if task_events_file does not exist, make it empty array [], 'x' would
% work too given how I wrote the scripts.
if ~exist('task_events_file', 'var') || ~ischar(task_events_file) || isempty(task_events_file)
    task_events_file=[]; 
elseif ~isfile(task_events_file)
    fprintf(['Task events file not found at ', task_events_file, '\n'])
    return
end

if ~exist('mel_fol', 'var') || ~ischar(mel_fol) || isempty(mel_fol)
    mel_fol=[]; % will end up using default
elseif ~isfolder(mel_fol)
    fprintf(['MELODIC folder not found at ', mel_fol, '\n'])
    mel_fol = [cicada_dir, '/sub-', subj_id, '/ses-', ses_id, '/', task_name, '/melodic'];
    fprintf(['Will run new melodic in new default folder: ' mel_fol, '\n'])
end

% First do all the checks to make sure everything seems legit!
% check that expected fmriprep directories exist first
fmriprep_func_dir = [fmriprep_dir, '/sub-', subj_id, '/ses-', ses_id, '/func'];
fmriprep_anat_dir = [fmriprep_dir, '/sub-', subj_id, '/ses-', anat_ses_id, '/anat'];
if ~isfolder(fmriprep_func_dir)
    fprintf(['Cannot find fmriprep func dir at ', fmriprep_func_dir, '\n'])
    return
elseif ~isfolder(fmriprep_anat_dir)
    fprintf(['Cannot find fmriprep anat dir at ', fmriprep_anat_dir, '\n'])
    return
else
    % first grab things from func folder
    cd(fmriprep_func_dir)
    funcfile_info = dir(['*', subj_id, '*', ses_id, '*', task_name, '*', 'space-MNI*preproc_bold.nii.gz']);
    funcfile = [funcfile_info.folder, '/', funcfile_info.name];
    funcmask_info = dir(['*', subj_id, '*', ses_id, '*', task_name, '*', 'space-MNI*brain_mask.nii.gz']);
    funcmask = [funcmask_info.folder, '/', funcmask_info.name];
    confounds_info =  dir(['*', subj_id, '*', ses_id, '*', task_name, '*', 'desc-confounds_timeseries.tsv']);
    confoundsfile = [confounds_info.folder, '/', confounds_info.name];

    % and now grab from best anat folder
    % Note, if you use a different structural other than T1w, you may need
    % to change that below
    cd(fmriprep_anat_dir)
    anatfile_info = dir(['*', subj_id, '*', anat_ses_id, '*', 'space-MNI*preproc_T1w.nii.gz']);
    anatfile = [anatfile_info.folder, '/', anatfile_info.name];
    anatmask_info = dir(['*', subj_id, '*', anat_ses_id, '*', 'space-MNI*desc-brain_mask.nii.gz']);
    anatmask = [anatmask_info.folder, '/', anatmask_info.name];
    gm_prob_info = dir(['*', subj_id, '*', anat_ses_id, '*', 'space-MNI*label-GM_probseg.nii.gz']);
    gm_prob = [gm_prob_info.folder, '/', gm_prob_info.name];
    wm_prob_info = dir(['*', subj_id, '*', anat_ses_id, '*', 'space-MNI*label-WM_probseg.nii.gz']);
    wm_prob = [wm_prob_info.folder, '/', wm_prob_info.name];
    csf_prob_info = dir(['*', subj_id, '*', anat_ses_id, '*', 'space-MNI*label-CSF_probseg.nii.gz']);
    csf_prob = [csf_prob_info.folder, '/', csf_prob_info.name];
end

output_dir = [cicada_dir, '/sub-', subj_id, '/ses-', ses_id, '/', task_name];

% Output everything to check:
fprintf(['\nfmriprep_dir: ', fmriprep_dir, '\n'])
fprintf(['cicada_dir: ', cicada_dir, '\n'])
fprintf(['subj_id: ', subj_id, '\n'])
fprintf(['ses_id: ', ses_id, '\n'])
fprintf(['task_name: ', task_name, '\n'])
fprintf(['ses_id for anat: ', anat_ses_id, '\n'])
fprintf(['task_events_file: ', task_events_file, '\n'])
fprintf(['compare_file: ', compare_file, '\n'])
fprintf(['melodic folder: ', mel_fol, '\n\n'])


% For melodic folder, we can just check if the melodic folder exists, does
% it have one of the final inputs?

Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, compare_file, task_events_file, mel_fol, anatfile, anatmask, gm_prob, wm_prob, csf_prob)

end