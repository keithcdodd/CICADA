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
smoothing_kernel = 0; % -1 would force default, which is 3mm smoothing kernel. Manuscript does no smoothing.
% More common would be something like 6 (i.e., 6 mm), but you can also just smooth the images yourself after CICADA if desired since smoothing is the last step!
fpass = []; % bandpassing in Hz. A great choice might be [0.008, 0.15]. Default is no bandpassing (like in the manuscript). [0.008, 0.1] is pretty commonplace.
detrended_degree = 2; % what degree of polynomial to detrend to. 2nd polynomial may be good as a bare minimum.
redo_melodic = 0;
sub_ids = {'132', '135', '140'};
ses_ids = {'01', '01', '01'}; % needs to be same length as sub_ids, since they will match up one for one
excludes = {'0', '0', '0'}; % if '1', you are excluding data from the dataset due to some error (e.g., scanner stopped halfway, incomplete image), see cicada_group_qc.m for more detail
outliers = {'0', '0', '0'}; % if '1', you are excluding data you have marked as an outlier (e.g., scanner issue that clearly corrupted image), see cicada_group_qc.m for more detail
adjusteds = {'0', '0', '0'}; % if '1', there are image(s) that you have performed Manual CICADA on (e.g., for some reason Auto CICADA was not very accurate for subject_01, so you performed manual ICA denoising through Manual CICADA and want to use that image instead.
task_event_files = {}; % the paths to the task event files, if this was task based.

cicada_group_qc(cicada_dir, group_qc_home, task_name, output_dirname, file_tag, smoothing_kernel, fpass, detrended_degree, redo_melodic, sub_ids, ses_ids, excludes, outliers, adjusteds, task_event_files)
