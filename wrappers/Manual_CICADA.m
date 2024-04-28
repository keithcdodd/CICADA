function Manual_CICADA(output_dir, compare_file, mel_fol, IC_manual_checker)
% Make this into a function that does name=value, and specify what inputs
% are required. 
% Necessary inputs: output_dir (e.g., cicada/subj-id/ses-id/task-id), funcfile, funcmask, and confoundsfile.
% Strongly suggested inputs: subject specific anatfile, anatmask, gm_prob,
% wm_prob, and csf_prob, as well as task_events_file (if task data).
% can skip variables with [] input


% need to make sure CICADA folder and subfolders are added to path!
Manual_CICADA_dir = fileparts(mfilename('fullpath')); % this gives current script path
cd([Manual_CICADA_dir, '/..'])
cicada_script_path = pwd;
cd(Manual_CICADA_dir);
addpath(genpath(cicada_script_path)); % add the basescripts to path if not already done 


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

if ~isfolder(output_dir)
    fprintf(['ERROR: Cannot find output_dir (task_dir) at ', output_dir, '\n'])
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

compare_file_2 = '';
% Now check for non-necessary variables
if ~exist('compare_file', 'var') || ~ischar(compare_file)
    fprintf('Will compare to standard 8 parameter and auto CICADA \n')
    compare_file='x'; % need it to be something that you can catch. 'x' is good.

    % It is good to compare manual to auto as well!
    compare_file_2_info = dir([output_dir, '/cleaned/*CICADA*auto*.nii.gz']);
    compare_file_2 = [compare_file_2_info.folder, '/', compare_file_2_info.name];
end

if ~exist('mel_fol', 'var') || ~ischar(mel_fol) || ~isfolder(mel_fol)
    fprintf('Will put Melodic in standard location \n')
    mel_fol = [output_dir, '/melodic'];
end

if ~isfolder(mel_fol)
    fprintf(['Cannot find melodic folder at ', mel_fol, '/n'])
    return
end

fprintf('\n\n')
%%%%%%%%%%%%%% ACTUALLY DO THE THINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('Running: CICADA_2_ManualLabeling \n')
cleaned_file = CICADA_2_ManualLabeling(output_dir, IC_manual_checker, mel_fol);
fprintf('Done with CICADA_2\n\n')

fprintf('Running: CICADA_3_QC \n')
CICADA_3_QC(cleaned_file, compare_file)
fprintf('Done with CICADA_3\n\n')

% And then if we are comparing to auto CICADA and the file is found,
% compare those too!
if ~isempty(compare_file_2) && isfile(compare_file_2)
    fprintf('Running: CICADA_3_QC to compare manual to auto \n')
    CICADA_3_QC(cleaned_file, compare_file_2)
    fprintf('Done with CICADA_3\n\n')
    % close figures
    close all
end

