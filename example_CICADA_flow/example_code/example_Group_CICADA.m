% After running example_Auto_CICADA, you can try running Group CICADA:

% Note: You will have to make sure the base directory is the example_CICADA
% directory, and that the path to the CICADA scripts is correct (and/or
% added)

CICADA_script_path = '/Users/keithdodd/GitHub/CICADA';
addpath(genpath(CICADA_script_path))
base_dir = '/Users/keithdodd/test_CICADA';
cicada_dir = [base_dir, '/bids_data/derivatives/cicada'];
group_qc_home = [base_dir, '/bids_data/derivatives/cicada_group_qc']; % general output directory
task_name = 'rest';
output_dirname = 'CICADA_auto_ses_01'; % unique name for the outputs specific to this group
file_tag = '_auto_'; % '_auto_' for Automatic CICADA, '_manual_' for Manual CICADA, because these tags are unique to auto and manual cicada files
voxelwise_scale = 0; % if set to 1, will perform voxelwise scaling after detrending, bandpassing, and smoothing
smoothing_kernel = 0; % -1 would force default, which is 3mm smoothing kernel (after detrending and bandpass filtering). Manuscript does no smoothing.
% More common would be something like 6 (i.e., 6 mm), but you can also just smooth the images yourself after CICADA if desired since smoothing is the last step!
fpass = []; % bandpassing in Hz, performed after detrending. A great choice might be [0.008, 0.15]. Default is no bandpassing (like in the manuscript). [0.008, 0.1] is pretty commonplace.
detrended_degree = 2; % what degree of polynomial to detrend to. 2nd polynomial may be good as a bare minimum.
redo_melodic = 0;
sub_ids = {'132', '135', '140'};
ses_id = {'01'}; % session id you are comparing across
exclude = {''}; % list of sub_ids you are excluding data from the dataset entirely (do not put into group CICADA folder) due to some error (e.g., scanner stopped halfway, incomplete image), see cicada_group_qc.m for more detail
outlier = {''}; % list of sub_ids you are excluding data you have marked as an outlier (e.g., scanner issue that clearly corrupted image), see cicada_group_qc.m for more detail
adjusted = {''}; % list of sub-ids that you have performed Manual CICADA on (e.g., for some reason Auto CICADA was not very accurate for subject_01, so you performed manual ICA denoising through Manual CICADA and want to use that image instead.
task_event_file = {''}; % path to the task event file relevant for all data (assuming same task, same timing).
% NOTE: if you have different task event files per subject/session, then
% you will need to edit task_event_files variable below



%%%%%%%%%%%%%% Should not need to change anything below, unless you have
%%%%%%%%%%%%%% different task event files per subject/session, then edit
%%%%%%%%%%%%%% task_event_files variable accordingly:%%%%%%%%%%%%%%%%%%%%%
% Change formating of variables appropriately
% Generate ses_ids
ses_ids = repmat(ses_id, size(sub_ids));

% Generate task event file list
if isempty(testing{1})
    task_event_files = {}; % proper formatting for CICADA
else
    % this is the variable you would edit if you have different task event
    % files per subject/session
    task_event_files = repmat(task_event_file, size(sub_ids));
end

% Initialize other output lists
excludes = cell(size(sub_ids));
outliers = cell(size(sub_ids));
adjusteds = cell(size(sub_ids));

% Fill in 0 or 1
for i = 1:length(sub_ids)
    excludes{i}  = num2str(double(ismember(sub_ids{i}, exclude)));
    outliers{i}  = num2str(double(ismember(sub_ids{i}, outlier)));
    adjusteds{i} = num2str(double(ismember(sub_ids{i}, adjusted)));
end


% Example formatting for three sub_ids in resting-state (no task):
% ses_ids = {'01', '01', '01'}; % needs to be same length as sub_ids, since they will match up one for one
% excludes = {'0', '0', '0'}; % if '1', you are excluding data from the dataset entirely (do not put into group CICADA folder) due to some error (e.g., scanner stopped halfway, incomplete image), see cicada_group_qc.m for more detail
% outliers = {'0', '0', '0'}; % if '1', you are excluding data you have marked as an outlier (e.g., scanner issue that clearly corrupted image), see cicada_group_qc.m for more detail
% adjusteds = {'0', '0', '0'}; % if '1', there are image(s) that you have performed Manual CICADA on (e.g., for some reason Auto CICADA was not very accurate for subject_01, so you performed manual ICA denoising through Manual CICADA and want to use that image instead.
% task_event_files = {}; % the paths to the task event files, if this was task based.


% Actually run Group CICADA!
cicada_group_qc(cicada_dir, group_qc_home, task_name, output_dirname, file_tag, voxelwise_scale, smoothing_kernel, fpass, detrended_degree, redo_melodic, sub_ids, ses_ids, excludes, outliers, adjusteds, task_event_files)
