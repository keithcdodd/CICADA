% Auto_CICADA.m should already set this all up for you as is, but in case you run into issues, here is the general code in order to run CICADA with fsl.

% This will almost certainly require modifications specific to your system, especially CICADA_software_dir and potentially your fsl path


CICADA_software_dir = '/home/GitHub/CICADA'; % edit to your actual path
addpath(genpath(CICADA_software_dir))
savepath


% Make sure fsl is set up correctly as well!
if (~contains(path, 'fsl*/etc/matlab')) || (~strcmp(getenv('FSLOUTPUTTYPE'), 'NIFTI_GZ'))
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

savepath