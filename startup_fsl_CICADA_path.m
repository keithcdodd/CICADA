% modify as needed for your CICADA, FSL, Matlab set up!
% you can then either run this yourself, or Auto_CICADA will run it for you too

% This will almost certainly require modifications specific to your system, especially CICADA_software_dir and potentially your fsl path

% NOTE: likely need to edit at least these next two lines
CICADA_software_dir = '/home/GitHub/CICADA'; % NOTE: edit to your actual path to CICADA
fslpath = '/usr/local/fsl'; % NOTE: edit to your actual path to fsl (which fsl)


addpath(genpath(CICADA_software_dir))

% Make sure fsl is set up correctly as well!
if isempty(regexp(path, '.*fsl.*/etc/matlab', 'once')) || (~strcmp(getenv('FSLOUTPUTTYPE'), 'NIFTI_GZ'))
    % FSL was not set up correctly, try to do that here
    fprintf('FSL with Matlab was not set up properly? Trying to do that for you now...\n')
    setenv( 'FSLDIR', fslpath ); % 
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

savepath