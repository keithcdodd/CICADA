function Auto_CICADA(output_dir, funcfile, funcmask, confoundsfile, compare_file, task_events_file, mel_fol, anatfile, anatmask, gm_prob, wm_prob, csf_prob)
% Make this into a function that does name=value, and specify what inputs
% are required. 
% Necessary inputs: output_dir (e.g., cicada/subj-id/ses-id/task-id), funcfile, funcmask, and confoundsfile.
% Strongly suggested inputs: subject specific anatfile, anatmask, gm_prob,
% wm_prob, and csf_prob, as well as task_events_file (if task data).
% can skip variables with [] input


% need to make sure CICADA folder and subfolders are added to path!
Auto_CICADA_dir = fileparts(mfilename('fullpath')); % this gives current script path
basescript_path = [Auto_CICADA_dir, '/../basescripts'];
cd(basescript_path);
basescript_dir=pwd;
cd(Auto_CICADA_dir);
addpath(basescript_dir); % add the basescripts to path if not already done 


% Make sure fsl is set up correctly as well!
if (~contains(path, 'fsl/etc/matlab')) || (~strcmp(getenv('FSLOUTPUTTYPE'), 'NIFTI_GZ'))
    % FSL was not set up correctly, try to do that here
    fprintf('FSL with Matlab was not set up properly? Trying to do that for you now...\n')
    setenv( 'FSLDIR', '/usr/local/fsl' );
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
    fsldir = getenv('FSLDIR');
    fsldirmpath = sprintf('%s/etc/matlab',fsldir);
    path(path, fsldirmpath);
    clear fsldir fsldirmpath;
end

% system path must also have fsl in it, otherwise will not work
curr_system_path = getenv('PATH');
fsldir = getenv('FSLDIR');
if ~contains(curr_system_path, fsldir)
    fprintf('Matlab System Path does not have fsl. Trying to add that for you now...\n') 
    new_system_path = [curr_system_path, ':', fsldir, '/bin'];
    setenv('PATH', new_system_path)
end
clear fsldir curr_system_path


% Now, check for required variables 
if ~exist('output_dir', 'var') || ~ischar(output_dir)
    fprintf('ERROR: Missing an output_dir specification or is not a character array!\n')
    return;
end

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
if ~exist('compare_file', 'var') || ~ischar(compare_file)
    fprintf('Will compare to standard 8 parameter \n')
    compare_file='x'; % need it to be something that you can catch. 'x' is good.
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
    fprintf('Will put Melodic in standard location \n')
    mel_fol='x'; % initialize it the same way the bash script would
end

fprintf('\n\n')
%%%%%%%%%%%%%% ACTUALLY DO THE THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CICADA_1_command = [basescript_dir, '/CICADA_1_MasksandICAs.sh', ' -o ', output_dir, ...
    ' -F ', funcfile, ' -f ', funcmask, ' -C ', confoundsfile, ...
    ' -A ', anatfile, ' -a ', anatmask, ' -g ', gm_prob, ...
    ' -w ', wm_prob, ' -c ', csf_prob, ' -m ', mel_fol, '< /dev/null'];

fprintf(['Running: ', CICADA_1_command, '\n'])
[status, cmdout_CICADA_1] = system(CICADA_1_command, '-echo');
fprintf('Done with CICADA_1\n\n')

fprintf('Running: CICADA_2_AutoLabeling \n')
tolerance = 4;
CICADA_2_AutoLabeling(output_dir, task_events_file, tolerance)
fprintf('Done with CICADA_2\n\n')

fprintf('Running: CICADA_3_QC \n')
cleaned_dir = [output_dir, '/cleaned_auto'];
CICADA_3_QC(cleaned_dir, compare_file)
fprintf('Done with CICADA_3\n\n')

end

