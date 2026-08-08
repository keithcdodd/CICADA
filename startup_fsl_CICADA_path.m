% CICADA / FSL MATLAB startup
%
% This script configures MATLAB to use CICADA and FSL.
%
% Auto_CICADA runs this script automatically. You may also run it yourself
% if you want to configure MATLAB for CICADA before starting an analysis.
%
%
% USER SETUP
% ==========
%
% In most cases, YOU DO NOT NEED TO EDIT THIS FILE.
%
% CICADA:
%   CICADA automatically determines its installation directory from the
%   location of this startup file. You do NOT need to enter the path to
%   your CICADA folder.
%
% FSL:
%   CICADA tries to locate FSL in the following order:
%
%     1. The optional manual FSL path specified below.
%     2. An already configured FSLDIR environment variable.
%     3. Common FSL installation locations.
%
%   Therefore, if FSL is already configured by your MATLAB startup file,
%   shell environment, or another setup script, CICADA will respect that
%   configuration rather than overwrite it.
%
%   If CICADA cannot locate FSL automatically, specify the full path to
%   your FSL installation below.
%
%   Example on macOS:
%       fslpath_override = '/Users/yourname/fsl';
%
%   Example on a server:
%       fslpath_override = '/path/to/software/fsl';
%
%   The specified directory should contain:
%       etc/fslconf/fsl.sh
%
%   Normally, leave the setting empty:
%       fslpath_override = '';
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OPTIONAL USER SETTING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fslpath_override = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normally no changes are needed below this line.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Locate CICADA automatically

% Because this startup file lives in the root CICADA directory, its own
% location provides the CICADA installation path.
CICADA_software_dir = fileparts(mfilename('fullpath'));

% Make the current CICADA checkout the preferred version on the MATLAB path.
addpath(genpath(CICADA_software_dir), '-begin');


%% Locate FSL

% Priority:
%   1. Explicit user override above
%   2. Existing valid FSLDIR
%   3. Common installation locations

if ~isempty(fslpath_override)

    % User explicitly supplied an FSL installation.
    fslpath = fslpath_override;

    if ~isfile(fullfile(fslpath, 'etc', 'fslconf', 'fsl.sh'))
        error('CICADA:FSLOverrideInvalid', ...
            ['The manually specified FSL directory does not appear to ' ...
             'contain a valid FSL installation:\n%s\n\n' ...
             'Expected to find:\n%s'], ...
            fslpath, ...
            fullfile(fslpath, 'etc', 'fslconf', 'fsl.sh'));
    end

else

    % Respect an existing FSLDIR if it points to a valid installation.
    fslpath = getenv('FSLDIR');

    fsl_valid = ...
        ~isempty(fslpath) && ...
        isfile(fullfile(fslpath, 'etc', 'fslconf', 'fsl.sh'));

    if ~fsl_valid

        % Try common installation locations.
        candidate_fslpaths = { ...
            fullfile(getenv('HOME'), 'fsl'), ...
            '/usr/local/fsl'};

        fslpath = '';

        for i = 1:numel(candidate_fslpaths)

            candidate = candidate_fslpaths{i};

            if isfile(fullfile(candidate, ...
                    'etc', 'fslconf', 'fsl.sh'))

                fslpath = candidate;
                break

            end

        end

        if isempty(fslpath)

            error('CICADA:FSLNotFound', ...
                ['CICADA could not locate a valid FSL installation.\n\n' ...
                 'Either:\n' ...
                 '  1. Configure FSLDIR before running CICADA, or\n' ...
                 '  2. Set fslpath_override near the top of\n' ...
                 '     startup_fsl_CICADA_path.m.\n\n' ...
                 'The FSL directory should contain etc/fslconf/fsl.sh.']);

        end

    end

end

setenv('FSLDIR', fslpath);


%% Configure FSL output format

if ~strcmp(getenv('FSLOUTPUTTYPE'), 'NIFTI_GZ')
    setenv('FSLOUTPUTTYPE', 'NIFTI_GZ');
end


%% Add FSL MATLAB support

fsldirmpath = fullfile(fslpath, 'etc', 'matlab');

if ~isfolder(fsldirmpath)

    error('CICADA:FSLMatlabMissing', ...
        ['CICADA found FSL at:\n%s\n\n' ...
         'but could not find its MATLAB support directory:\n%s'], ...
        fslpath, fsldirmpath);

end

call_fsl_path = which('call_fsl');

% Add the MATLAB support belonging to the selected FSL installation.
% '-begin' ensures that an older/stale FSL MATLAB path does not take
% precedence over the installation selected above.
if isempty(call_fsl_path) || ...
        ~contains(call_fsl_path, fsldirmpath)

    addpath(fsldirmpath, '-begin');

end


%% Make FSL command-line programs available to CICADA shell scripts

% Standard FSL installations may expose executable entry points in
% FSLDIR/share/fsl/bin and/or FSLDIR/bin. Add whichever exist.

fsl_command_dirs = { ...
    fullfile(fslpath, 'share', 'fsl', 'bin'), ...
    fullfile(fslpath, 'bin')};

valid_fsl_command_dirs = {};

for i = 1:numel(fsl_command_dirs)

    if isfolder(fsl_command_dirs{i})
        valid_fsl_command_dirs{end+1} = fsl_command_dirs{i}; %#ok<SAGROW>
    end

end

if isempty(valid_fsl_command_dirs)

    error('CICADA:FSLBinMissing', ...
        ['CICADA found FSL at:\n%s\n\n' ...
         'but could not locate its command-line programs.'], ...
        fslpath);

end

curr_system_path = getenv('PATH');
path_entries = strsplit(curr_system_path, pathsep);

% Add FSL command locations to the beginning of MATLAB's system PATH.
% This is necessary because CICADA also calls FSL from shell scripts.
for i = numel(valid_fsl_command_dirs):-1:1

    this_dir = valid_fsl_command_dirs{i};

    if ~any(strcmp(path_entries, this_dir))

        curr_system_path = ...
            [this_dir, pathsep, curr_system_path];

        path_entries = [{this_dir}, path_entries];

    end

end

setenv('PATH', curr_system_path);


%% Report the resolved configuration

fprintf('CICADA directory: %s\n', CICADA_software_dir);
fprintf('FSL directory:    %s\n', fslpath);


%% Clean temporary startup variables

clear fslpath_override fslpath fsl_valid ...
    candidate_fslpaths candidate i ...
    fsldirmpath call_fsl_path ...
    fsl_command_dirs valid_fsl_command_dirs ...
    curr_system_path path_entries this_dir;