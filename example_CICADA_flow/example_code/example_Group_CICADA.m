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
voxelwise_scale_flag = 0; % if set to 1, will perform voxelwise scaling after detrending, bandpassing, and smoothing
smoothing_kernel = 0; % -1 would force default, which is 3mm smoothing kernel (after detrending and bandpass filtering). Manuscript does no smoothing.
% More common would be something like 6 (i.e., 6 mm), but you can also just smooth the images yourself after CICADA if desired since smoothing is the last step!

compare_tag = '8p'; % comparing to 8 parameter regression is standard, 
% but this just needs to be an identifiable tag of other cleaned data in your cleaned dir that you want to compare to. 
% e.g.,: '6p', '8p', '9p', '12p', '16p', '18p', '24p', '28p', '30p', '32p', '36p', or
% whatever fits what you already have in your cleaned_dirs!

fpass = []; % bandpassing in Hz, performed after detrending. A great choice might be [0.008, 0.15]. Default is no bandpassing (like in the manuscript). [0.008, 0.1] is pretty commonplace.
detrended_degree = 2; % what degree of polynomial to detrend to. 2nd polynomial may be good as a bare minimum.
redo_melodic = 0;
sub_ids = {'132', '135', '140'};
ses_id = {'01'}; % session id you are comparing across
exclude = {''}; % list of sub_ids you want to fully exclude, for example had clear scanning issues or some other reason why it should be excluded from consideration
adjusted = {''}; % list of sub-ids that you have performed Manual CICADA on (e.g., for some reason Auto CICADA was not very accurate for subject_01, so you performed manual ICA denoising through Manual CICADA and want to use that image instead.
task_event_file = {''}; % path to the task event file relevant for all data (assuming same task, same timing).
% NOTE: if you have different task event files per subject/session, then
% you will need to edit task_event_files variable below


% Optional:
% Subject level data table inclusion so that Group CICADA can include
% other information in final group_qc_table.csv such as outcome measures
% and covariates like age, sex, etc. Whatever you might want to include in
% future analyses! Just to make things easier later :)
subject_level_data_table = table(); % optional, If no table, just keep subject level data table as an empty table() 
subject_level_data_table.subject = char(sub_ids); % this column is required if you use this table. Should have at least all your sub_ids you are including above, can have more.
subject_level_data_table.age = ['40'; '20'; '33']; % you can include any and as many subject level data columns you want.
subject_level_data_table.sex = ['F'; 'M'; 'F'];
% alternatively, you could have made a .csv table already and just use
% readtable() function! This is often much easier especially with larger
% data sizes. The .csv table just has to at least have the sub_ids, and
% only one instance of each! Can have more sub_ids within it without issue
% (e.g., sub_ids you end up not measuring/using for example)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%% Should not need to change anything below this, unless you have
%%%%%%%%%%%%%% different task event files per subject/session, then edit
%%%%%%%%%%%%%% task_event_files variable below accordingly:%%%%%%%%%%%%%%%%%%%%%

% Change formating of variables appropriately
% Generate ses_ids
ses_ids = repmat(ses_id, size(sub_ids));

% Generate task event file list
if isempty(task_event_file{1})
    task_event_files = {}; % proper formatting for CICADA
else
    % this is the variable you would edit if you have different task event
    % files per subject/session
    task_event_files = repmat(task_event_file, size(sub_ids));
end

% Initialize other output lists
excludes = cell(size(sub_ids));
adjusteds = cell(size(sub_ids));

% Fill in 0 or 1
for i = 1:length(sub_ids)
    excludes{i}  = num2str(double(ismember(sub_ids{i}, exclude)));
    adjusteds{i} = num2str(double(ismember(sub_ids{i}, adjusted)));
end


% Example formatting for three sub_ids in resting-state (no task):
% ses_ids = {'01', '01', '01'}; % needs to be same length as sub_ids, since they will match up one for one
% excludes = {'0', '0', '0'}; % if '1', you are excluding data from the dataset entirely (do not put into group CICADA folder) due to some error (e.g., scanner stopped halfway, incomplete image), see cicada_group_qc.m for more detail
% adjusteds = {'0', '0', '0'}; % if '1', there are image(s) that you have performed Manual CICADA on (e.g., for some reason Auto CICADA was not very accurate for subject_01, so you performed manual ICA denoising through Manual CICADA and want to use that image instead.
% task_event_files = {}; % the paths to the task event files, if this was task based.


% Actually run Group CICADA!
cicada_group_qc(cicada_dir, group_qc_home, task_name, output_dirname, file_tag, voxelwise_scale_flag, smoothing_kernel, fpass, detrended_degree, redo_melodic, sub_ids, ses_ids, excludes, adjusteds, compare_tag, task_event_files, subject_level_data_table)
