function Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, redo_mel, mel_fol, compare_file, task_events_file, anatfile, anatmask, gm_prob, wm_prob, csf_prob, tolerance, despike)
% Make this into a function that does name=value, and specify what inputs
% are required.
% Necessary inputs: output_dir (e.g., cicada/subj-id/ses-id/task-id), funcfile, funcmask, and confoundsfile.
% Strongly suggested inputs: subject specific anatfile, anatmask, gm_prob,
% wm_prob, and csf_prob, as well as task_events_file (if task data).
% can skip variables with [] input
% compare_file allows you to compare the denoising to any file you want,
% just needs the full path. Alternatively, you can set compare_file to the
% characters of either '6p', '8p', '9p', '12p', '16p', '18p',
% '24p', '28p', '30p', '32p', or '36p' to use one of CICADA's
% standard confound regressions for comparison.
% despike 1 will do some light despiking of the data before diving into
% CICADA, despike 0 (or left blank) will not

% need to make sure CICADA folder and subfolders are added to path!
Auto_CICADA_dir = fileparts(mfilename('fullpath')); % this gives current script path
basescript_path = [Auto_CICADA_dir, '/../basescripts'];
cd(basescript_path);
basescript_dir=pwd;
cd(Auto_CICADA_dir);
addpath(basescript_dir); % add the basescripts to path if not already done

% Run the shared portable CICADA/FSL startup. It respects an existing
% valid FSLDIR and documents an optional manual FSL-path override.
if ~isfile([Auto_CICADA_dir, '/../startup_fsl_CICADA_path.m'])
	fprintf(['Cannot find ', Auto_CICADA_dir, '/../startup_fsl_CICADA_path.m !!!\n'])
	return
end

run([Auto_CICADA_dir, '/../startup_fsl_CICADA_path.m']) % sets up FSL and Matlab, if not done already

% Now, check for required variables
if ~exist('output_dir', 'var') || ~ischar(output_dir)
    fprintf('ERROR: Missing an output_dir specification or is not a character array!\n')
    return;
end

if ~exist('despike', 'var') || ~ismember(despike, [0 1])
    despike = 0;
end

% now create a diary log in output dir and record everything per output:
log_dir = [output_dir, '/logs'];
if isfolder(log_dir)
    rmdir(log_dir, 's')
end
mkdir(log_dir)
date_label = char(datetime('today'));
diaryfile = [log_dir, '/Auto_CICADA_', date_label, '_log'];
diary(diaryfile)

if ~exist('funcfile', 'var') || ~ischar(funcfile)
    fprintf('ERROR: Missing a funcfile specification or is not a character array!\n')
    return;
end

if ~exist('funcmask', 'var') || ~ischar(funcmask)
    fprintf('ERROR: Missing a funcmask specification or is not a character array!\n')
    return;
end

if ~exist('confoundsfile', 'var') || ~ischar(confoundsfile)
    fprintf('ERROR: Missing a confoundsfile or is not a character array!\n')
    return;
end

% Now check for non-necessary variables
if ~exist('redo_mel', 'var') || (redo_mel ~= 0 && redo_mel ~=1)
    redo_mel = 0; % default is to not redo melodic
end

compare_tag = '8p'; % set a default just in case
if ~exist('compare_file', 'var') || ~ischar(compare_file) || isempty(compare_file)
    compare_file = '';
elseif ~isfile(compare_file)
    % it is a char array, but not a valid file
    % check to see if it is a valid tag instead!
    valid_tags = {'6p', '8p', '9p', '12p', '16p', '18p', '24p', '28p', '30p', '32p', '36p'}; % for compare file if you want to use an inbuilt one!
    if ~ismember(compare_file, valid_tags)
        fprintf('Not a valid tag or file... Will compare to 8p regression.\n')
        compare_file = '';
    else
        fprintf(['Will compare to ', compare_file, ' regression!\n'])
    end
end

% at this point, compare_file is either a valid file, or a valid
% compare_tag, or empty char array

if ~isfile(compare_file)
    compare_tag = compare_file;

    if isempty(compare_tag)
        compare_tag = '8p'; % to the default! Because it was not a valid file or tag
    end
end

if ~exist('task_events_file', 'var')
    task_events_file='x';
end

if ~exist('anatfile', 'var') || ~ischar(anatfile)
    anatfile='x'; % initialize it the same way the bash script would
end

if ~exist('anatmask', 'var') || ~ischar(anatmask)
    anatmask='x'; % initialize it the same way the bash script would
end

if ~exist('gm_prob', 'var') || ~ischar(gm_prob)
    gm_prob='x'; % initialize it the same way the bash script would
end

if ~exist('wm_prob', 'var') || ~ischar(wm_prob)
    wm_prob='x'; % initialize it the same way the bash script would
end

if ~exist('csf_prob', 'var') || ~ischar(csf_prob)
    csf_prob='x'; % initialize it the same way the bash script would
end

if ~exist('mel_fol', 'var') || ~ischar(mel_fol)
    fprintf('Melodic will be in the standard location \n')
    mel_fol = [output_dir, '/melodic']; % give it default location
end

% change mel_fol if redo melodic is a 1
if redo_mel == 1
    fprintf('Will run a new Melodic instance.\n')
    mel_fol = 'x'; % initilaize it the same was as the basescript would
end

% if applicable, despike first
if despike == 1

    fprintf('Lightly Despiking the Data...\n')

    % If no GM probability map was supplied, use the same default MNI GM
    % probability map that CICADA_1 would normally select for gm_prob='x'.
    if strcmp(gm_prob, 'x')
        gm_prob = [Auto_CICADA_dir, ...
            '/../templates/mni_icbm152_nlin_asym_09c/', ...
            'mni_icbm152_gm_tal_nlin_asym_09c.nii.gz'];

        fprintf('Using CICADA default MNI GM probability map for despiking.\n')
    end

    if ~isfile(gm_prob)
        fprintf(['ERROR: Cannot find GM probability map for despiking: ', ...
            gm_prob, '\n'])
        return
    end

    cut2 = 6; % CICADA-tested localedit spikiness threshold -> cut2 = 6.
          % Lower values edit more voxel-time samples; higher values edit fewer.

    [funcfile, ~, despike_qc] = despike_fMRI( ...
        funcfile, gm_prob, cut2, output_dir);

    fprintf(['Despiking complete: cut2=%.6g; %d voxel-time samples edited ', ...
        '(%.6f%% of GM voxel-time).\n'], ...
        despike_qc.cut2, ...
        despike_qc.n_altered_samples, ...
        100 * despike_qc.fraction_gm_voxel_time_samples);

    if redo_mel == 0 && ischar(mel_fol) && isfolder(mel_fol)
        warning(['Despiking is enabled while an existing MELODIC folder may ', ...
            'be reused. Confirm that the existing MELODIC decomposition was ', ...
            'generated from the same despiked functional data. Otherwise ', ...
            'rerun MELODIC.']);
    end
end

% Keep a diary log of output from here to start CICADA

fprintf('\n\n')
%%%%%%%%%%%%%% ACTUALLY DO THE THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CICADA_1_command = [basescript_dir, '/CICADA_1_MasksandICAs.sh', ' -o ', output_dir, ...
    ' -F ', funcfile, ' -f ', funcmask, ' -C ', confoundsfile, ...
    ' -A ', anatfile, ' -a ', anatmask, ' -g ', gm_prob, ...
    ' -w ', wm_prob, ' -c ', csf_prob, ' -m ', mel_fol, '< /dev/null'];

fprintf(['Running: ', CICADA_1_command, '\n'])
[status, cmdout_CICADA_1] = system(CICADA_1_command, '-echo');
if status ~= 0
    fprintf(2, ['ERROR: CICADA_1 failed with status ', ...
        num2str(status), '. Stopping this CICADA run.\n'])

    if ~isempty(strtrim(cmdout_CICADA_1))
        fprintf(2, '\nCICADA_1 command output:\n%s\n', cmdout_CICADA_1);
    end
    diary off
    return
end

fprintf('Done with CICADA_1\n\n')

fprintf('Running: CICADA_2_AutoLabeling \n')
if ~exist('tolerance', 'var')
    tolerance = 5;
end
if isfile(compare_file)
    cleaned_file = CICADA_2_AutoLabeling(output_dir, task_events_file, [], tolerance, mel_fol); % allow one to use another inbuilt compare file if desired, will still run 8p denoising regardless
else
    cleaned_file = CICADA_2_AutoLabeling(output_dir, task_events_file, compare_tag, tolerance, mel_fol); % allow one to select to compare to any of the inbuilt comparisons if desired
end
fprintf('Done with CICADA_2\n\n')

fprintf('Running: CICADA_3_QC \n')
if isfile(compare_file)
    CICADA_3_QC(cleaned_file, compare_file) % compare to a given file, which could be an inbuilt one
else
    % could be a valid tag, or empty char array
    CICADA_3_QC(cleaned_file, compare_tag) % compare to tag or default 8p
end

fprintf('Done with CICADA_3\n\n')
% close figures from CICADA_3_QC
close all

diary off

end

