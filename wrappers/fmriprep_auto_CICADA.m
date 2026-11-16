function fmriprep_auto_CICADA(fmriprep_dir, cicada_dir, sub_id, ses_id, task_name, anat_ses_id, redo_mel, mel_fol, ...
    task_events_file, compare_file, tolerance, despike, task_timing_options, template_space, template_res)
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

% Optional spatial template.
if nargin < 14 || isempty(template_space)
    template_space = 'MNI152NLin6Asym';
end

% Optional functional resolution label.
% Empty means automatically accept a unique matching resolution.
% Examples: '2', '02', '1', 'native'.
if nargin < 15 || isempty(template_res)
    template_res = '';
end

% Anatomical and functional files must use the same template,
% but need not have identical spatial sampling resolutions.
space_pattern = ['space-' template_space '_*'];

if isempty(template_res)
    func_space_pattern = space_pattern;
else
    template_res = char(string(template_res));
    func_space_pattern = ...
        ['space-' template_space '_res-' template_res '_*'];
end

sub_dir = fullfile(fmriprep_dir, ['sub-' sub_id]);

if ~isempty(ses_id)
    session_dir = fullfile(sub_dir, ['ses-' ses_id]);

    assert(isfolder(session_dir), ...
        'Requested session directory not found: %s', session_dir);

    has_ses = true;
else
    has_ses = false;
    anat_ses_id = '';
end

% Now check for non-necessary variables
if ~exist('redo_mel', 'var') || (redo_mel ~= 0 && redo_mel ~=1)
    redo_mel = 0; % default is to not redo melodic
end

if ~exist('despike', 'var') || (despike ~= 0 && despike ~=1)
    despike = 0; % default is to not despike
end

% Optional timing options; old calls remain valid.
if nargin < 13 || isempty(task_timing_options)
    task_timing_options = struct();
end

assert(isstruct(task_timing_options) && isscalar(task_timing_options), ...
    'CICADA:TimingOptions', ...
    'task_timing_options must be a scalar struct.');

if ~exist('compare_file', 'var') || ~ischar(compare_file) || isempty(compare_file)
    compare_file = '';
    compare_file_record = 'Standard 8 parameter compare';
elseif ~isfile(compare_file)
    % check to see if it is a valid tag instead!
    valid_tags = {'6p', '8p', '9p', '12p', '16p', '18p', ...
    '24p', '28p', '30p', '32p', '36p'};
    if ~ismember(compare_file, valid_tags)
        fprintf('Not a valid tag or file... Will compare to 8p regression.\n')
        compare_file = '';
        compare_file_record = 'Standard 8 parameter compare';
    else
        fprintf(['Will compare to ', compare_file, ' regression!\n'])
        compare_file_record = [compare_file, ' compare'];
    end
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

% Anatomical derivatives may be stored at subject level.
if ~isfolder(fmriprep_anat_dir)
    fallback_anat_dir = fullfile(fmriprep_dir, ...
        ['sub-' sub_id], 'anat');

    if isfolder(fallback_anat_dir)
        fmriprep_anat_dir = fallback_anat_dir;
    end
end

if ~isfolder(fmriprep_func_dir)
    fprintf('Cannot find fMRIPrep functional directory: %s\n', ...
        fmriprep_func_dir);
    return;
else
    % first grab things from func folder
    cd(fmriprep_func_dir)

    % Construct the exact subject/session/task prefix.
    if has_ses
        run_prefix = ['sub-' sub_id '_ses-' ses_id ...
            '_task-' task_name];
    else
        run_prefix = ['sub-' sub_id '_task-' task_name];
    end

    % Match the selected template and optional resolution.
    funcfile_info = dir(fullfile(fmriprep_func_dir, ...
        [run_prefix '_*' func_space_pattern ...
        'desc-preproc_bold.nii.gz']));

    assert(numel(funcfile_info) == 1, ...
        ['Expected exactly one preprocessed BOLD file. ' ...
         'Check task_name, mni_space, and mni_res.']);

    funcfile = fullfile(funcfile_info.folder, ...
        funcfile_info.name);

    % Derive the corresponding mask from the exact BOLD filename.
    funcmask = strrep(funcfile, ...
        'desc-preproc_bold.nii.gz', ...
        'desc-brain_mask.nii.gz');

    assert(isfile(funcmask), ...
        'Missing functional brain mask: %s', funcmask);

    % Confounds must correspond to this acquisition.
    bold_name = string(funcfile_info.name);
    assert(contains(bold_name, "_space-"), ...
        'Cannot identify BIDS space entity.');

    confounds_prefix = extractBefore(bold_name, "_space-");

    confounds_info = dir(fullfile(fmriprep_func_dir, ...
        char(confounds_prefix + "_desc-confounds_*.tsv")));

    assert(numel(confounds_info) == 1, ...
        'Expected exactly one matching confounds TSV.');

    confoundsfile = fullfile(confounds_info.folder, ...
        confounds_info.name);

    if ~isfile(confoundsfile)
        fprintf(['Cannot find confounds_file at ', confoundsfile, '\n'])
        return;
    end

    % Search session-level anatomy first, then subject-level anatomy.
    anat_dirs = unique({ ...
        fmriprep_anat_dir, ...
        fullfile(fmriprep_dir, ['sub-' sub_id], 'anat')}, ...
        'stable');

    anatfile_info = [];

    for d = 1:numel(anat_dirs)

        this_dir = anat_dirs{d};

        if ~isfolder(this_dir)
            continue;
        end

        % Prefer anatomy associated with the requested session.
        if ~isempty(anat_ses_id)
            candidate = dir(fullfile(this_dir, ...
                ['sub-' sub_id '_ses-' anat_ses_id '_' ...
                space_pattern 'desc-preproc_T1w.nii.gz']));

            if ~isempty(candidate)
                assert(numel(candidate) == 1, ...
                    'Multiple matching session-level T1 images.');
                anatfile_info = candidate;
                break;
            end
        end

        % Otherwise look for subject-level T1 without ses- entity.
        candidate = dir(fullfile(this_dir, ...
            ['sub-' sub_id '_' space_pattern ...
            'desc-preproc_T1w.nii.gz']));

        if ~isempty(candidate)
            assert(numel(candidate) == 1, ...
                'Multiple matching subject-level T1 images.');
            anatfile_info = candidate;
            break;
        end
    end

    assert(numel(anatfile_info) == 1, ...
        'Could not uniquely identify T1 in the requested space.');

    fmriprep_anat_dir = anatfile_info.folder;

    anatfile = fullfile(fmriprep_anat_dir, ...
        anatfile_info.name);

    anat_prefix = strrep(anatfile_info.name, ...
        'desc-preproc_T1w.nii.gz', '');

    anatmask = fullfile(fmriprep_anat_dir, ...
        [anat_prefix 'desc-brain_mask.nii.gz']);

    gm_prob = fullfile(fmriprep_anat_dir, ...
        [anat_prefix 'label-GM_probseg.nii.gz']);

    wm_prob = fullfile(fmriprep_anat_dir, ...
        [anat_prefix 'label-WM_probseg.nii.gz']);

    csf_prob = fullfile(fmriprep_anat_dir, ...
        [anat_prefix 'label-CSF_probseg.nii.gz']);

    required_anat = {anatfile, anatmask, gm_prob, wm_prob, csf_prob};

    for k = 1:numel(required_anat)
        assert(isfile(required_anat{k}), ...
            'Missing anatomical derivative: %s', required_anat{k});
    end
end

if has_ses == 1
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-', ses_id, '/', task_name];
else
    % just put in a ses-01 folder for consistency
    output_dir = [cicada_dir, '/sub-', sub_id, '/ses-01/', task_name];
end

if ~exist('tolerance', 'var')
    tolerance = 5;
end

% Output everything to check:
fprintf(['\nfmriprep_dir: ', fmriprep_dir, '\n'])
fprintf(['cicada_dir: ', cicada_dir, '\n'])
fprintf(['sub_id: ', sub_id, '\n'])
fprintf(['ses_id: ', ses_id, '\n'])
fprintf(['task_name: ', task_name, '\n'])
fprintf(['anat folder location: ', fmriprep_anat_dir, '\n'])
fprintf(['task_events_file: ', task_events_file_record, '\n'])
fprintf(['compare_file: ', compare_file_record, '\n'])
fprintf(['melodic folder: ', mel_fol, '\n\n'])
fprintf(['despike: ', num2str(despike), '\n'])


% For melodic folder, we can just check if the melodic folder exists, does
% it have one of the final inputs?

Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, redo_mel, mel_fol, compare_file, task_events_file, anatfile, anatmask, gm_prob, wm_prob, csf_prob, tolerance, despike, task_timing_options)

end