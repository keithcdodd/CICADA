% modify as needed for your CICADA, FSL, Matlab set up!
% you can then either run this yourself, or Auto_CICADA will run it for you too

% This will almost certainly require modifications specific to your system, especially CICADA_software_dir and potentially your fsl path

% NOTE: likely need to edit at least these next two lines
CICADA_software_dir = '/Users/keithdodd/GitHub/CICADA'; % NOTE: edit to your actual path to CICADA
fslpath = '/usr/local/fsl'; % NOTE: edit to your actual path to fsl directory that your terminal uses which should have the /etc/matlab within it!

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% edited 03/10/26 to hopefully get less warnings regarding savepath if on
% server.
path_edited = 0;
% see if it is already added by checking for presence of a key function
if isempty(which('CICADA_2_AutoLabeling')) || ~contains(path, CICADA_software_dir)
    addpath(genpath(CICADA_software_dir));
    path_edited = 1;
end

% Make sure fsl is set up correctly as well!
if (~strcmp(getenv('FSLOUTPUTTYPE'), 'NIFTI_GZ'))
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
end

if ~strcmp(getenv('FSLDIR'), fslpath)
    setenv('FSLDIR', fslpath);
end

% Add FSL MATLAB support
fsldir = getenv('FSLDIR');
fsldirmpath = fullfile(getenv('FSLDIR'),'etc','matlab');
if isempty(regexp(path, '.*fsl.*/etc/matlab', 'once'))
    % FSL was not set up correctly, try to do that here
    fprintf('FSL with Matlab was not set up properly? Trying to do that for you now...\n')
    addpath(fsldirmpath,'-begin') % safer than path command, add at beginning of path
    path_edited = 1;
end

% system path must also have fsl in it, otherwise will not work
curr_system_path = getenv('PATH');
if ~contains(curr_system_path, fsldir)
    fprintf('Matlab System Path does not have fsl. Trying to add that for you now...\n') 
    new_system_path = [curr_system_path, pathsep, fullfile(fsldir,'bin')];
    setenv('PATH', new_system_path)
end
clear fsldir curr_system_path fsldirmpath;

if path_edited == 1
    savepath
end
