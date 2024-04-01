function [cleaned_file] = CICADA_2_ManualLabeling(output_dir, IC_manual_checker, mel_fol)
% 
% NOTE: NEED TO MAKE IC_manual_checker.csv BEFORE RUNNING THIS FUNCTION:
% To make IC_manual_checker.csv, open IC_auto_checker.csv in
% ic_auto_selection folder (assuming you have run the auto version already). 
% Use melodic report, IC_probabilities, and any
% other information you might want to make decisions to edit the Potential
% Signal Labels column in IC_auto_checker.csv (1 means signal, 0 is noise) 
% [we suggest using Griffanti et al. 2017 Hand classification of fMRI ICA 
% noise components paper as a reference]. and then save as 
% IC_manual_checker.csv in that same location (ic_auto_selection folder). 
% THEN you can run this function.
% IC_manual_checker is the IC_manual_checker.csv

% output_dir must be the same one as for the first script, so that this can
% find all the relevant files. Again, it is best if this is in the format
% of CICADA_dir/subj-id/ses-id/task-id (so, a task_dir, really)

% 
%
% IC_manual_checker (optional): IC_manual_checker.csv file - if not given, will look for it in ic_auto_selection
% folder

% check output_dir exists 
if ~isfolder(output_dir)
    fprintf('ERROR: output_dir does not exist!')
    return;
end

% read in task name, session name, and subject name based on default folder
% structure
cd(output_dir)
task_dir = output_dir;

if ~exist('mel_fol', 'var') || isempty(mel_fol) || ~isfolder(mel_fol)
    mel_fol = [task_dir, '/melodic'];
end

if ~isfolder(mel_fol)
    fprintf(['Cannot find melodic folder at ' mel_fol, '\n'])
    return;
end

if ~exist('IC_manual_checker', 'var')
    % If not given assume it is in its normal spot and look for it
    if isfile([output_dir, '/ic_auto_selection/IC_manual_checker.csv'])
        IC_manual_checker = [output_dir, '/ic_auto_selection/IC_manual_checker.csv'];
    elseif isfile([output_dir, '/ic_manual_selection/IC_manual_checker.csv'])
            IC_manual_checker = [output_dir, '/ic_manual_selection/IC_manual_checker.csv'];
    else
        fprintf('Cannot find IC_manual_checker. Exiting...\n')
        return
    end
else
    if ~isfile(IC_manual_checker)
        fprintf(['Could not find ', IC_manual_checker, ' ... Exiting...\n'])
        return
    end
end

% set up
manual_table = readtable(IC_manual_checker);
if ~isfile([output_dir, '/ic_auto_selection/DecisionVariables_Auto.mat'])
    fprintf('Cannot find ic_auto_selection/DecisionVariables_Auto.mat ... Exiting ...\n')
    return
end

load([output_dir, '/ic_auto_selection/DecisionVariables_Auto.mat'], 'output_dir', 'subject_id', 'session_id', 'task_id', 'Data', 'Tables', 'Results')

Results_Auto = Results; % Save auto results struct
% Create a new Results struct
Results = struct;

% Manual table should have all the same as auto table except for the
% SignalLabel column (which should be updated):
auto_table = Results_Auto.IC_checker_table;

% do some extra work to make sure we have the same columns as auto table
% First, make sure the sort of the manual table is the same as the auto
% table:
temp_manual_table = manual_table;
temp_auto_table = auto_table;
[B_manual, ~] = sortrows(temp_manual_table, 'PotentialICs', 'ascend'); % put in order from IC 1:end
[~, index_auto] = sortrows(temp_auto_table, 'PotentialICs', 'ascend'); % put in order from IC 1:end
fixed_manual_table = table(); % initialize as a table
fixed_manual_table(index_auto,:) = B_manual; % OK, now manual table PotentialICs order matches auto table even if it did not before
% Now, create a new IC Checker that is in the same form as auto, but
% matches decisions from manual
temp_table = auto_table; % start it as a match to auto_table
temp_table.SignalLabel = fixed_manual_table.SignalLabel; % update the signal label to manual
IC_checker_table = temp_table; % update name to IC_checker_table
% IC_checker_table should be the same as the imported manual table, but in
% case you changed other columns in the auto checker, this makes sure it
% has the same auto columns, just updates the SignalLabel column
signal_idx = IC_checker_table.PotentialICs; % This just gives us the ordered list of ICs

% input into Manual_Results
Results.IC_checker_table = IC_checker_table;

% Note what the changes are in Signal Label Column
IC_selection_changes = IC_checker_table(IC_checker_table.SignalLabel ~= auto_table.SignalLabel, :); % this will list the signal labels that manual has, that auto selection had the opposite on

writetable(IC_selection_changes, 'changed_signal_label_manual_ICs.csv')

% input into Manual_Results
Results.IC_selection_changes = IC_selection_changes;

% Now do the things, cut off the table to not include the classification
% strings
cut_off_table = IC_checker_table(:,{'PotentialICs', 'SignalLabel', 'HighSignalLabel', 'HighNoiseLabel'});
manual_array = table2array(cut_off_table);

% import some variables from Auto Results:
noise_indices = Results_Auto.noise_indices;
ICs = Results_Auto.ICs;
IC_exp_var = table2array(Tables.IC_exp_var_table);

% update noise_indices
noise_indices(manual_array(manual_array(:,2) == 0,1)) = 1;
noise_indices(manual_array(manual_array(:,2) == 1,1)) = 0;
noise_indices = logical(noise_indices);

% Create final noise and signal ICs and indices
noise_ICs = ICs(noise_indices);
signal_indices = logical(~noise_indices);
signal_ICs = ICs(signal_indices);

% Update manual results struct
Results.signal_indices = signal_indices;
Results.signal_ICs = signal_ICs;
Results.noise_indices = noise_indices;
Results.noise_ICs = noise_ICs;

% Make it easy to compare groups
eval_table = Results_Auto.feature_table_norm(int32(manual_array(:,1)), :);
Results.eval_table = eval_table;

% Compare before and after selection
feature_relative_array = table2array(Results_Auto.feature_relative_table(:,2:end));
before = sum(feature_relative_array .* IC_exp_var);
IC_exp_var_signal = IC_exp_var(signal_ICs);
IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
after = sum(feature_relative_array(signal_ICs,:) .* IC_exp_var_signal,1);
percent_change = 100.*(after ./ before - 1); % rounded to be easier to read
compare_cleaning = array2table([before; after; percent_change]', 'RowNames', Results_Auto.feature_relative_table.Properties.VariableNames(2:end), 'VariableNames', {'Before', 'After', 'Percent_Change'});
Results.compare_cleaning = compare_cleaning;

% estimate effective degrees of freedom (somewhat) from CADICA alone via proportions:
Results.percent_variance_kept = sum(IC_exp_var(signal_ICs)); % probably want > 0.1 (10%)
Results.num_ICs_kept = length(signal_ICs);
Results.num_ICs_total = length(ICs);
Results.IC_exp_var = IC_exp_var;
Results.ICs = ICs;

% Create images to review signal and noise components from
% probabilities ICs
all_prob = niftiread([mel_fol, '/ICprobabilities.nii.gz']);
gm_prob = niftiread([task_dir, '/region_masks/GM_prob.nii.gz']);
notgm_prob = niftiread([task_dir, '/region_masks/NotGM_prob.nii.gz']);
gm_bin = gm_prob > 0.75;
notgm_bin = notgm_prob > 0.75;
noise_prob = all_prob(:,:,:,noise_ICs);
noise_prob_info = niftiinfo([mel_fol, '/ICprobabilities.nii.gz']);
signal_prob = all_prob(:,:,:,signal_ICs);

% grab noise probabilities, update header, and write, same for
% signal
noise_prob_info.ImageSize = size(noise_prob);
noise_prob_info.DisplayIntensityRange = [0.9 1];
noise_prob_info.Datatype = 'single';

signal_prob_info = noise_prob_info;
signal_prob_info.ImageSize = size(signal_prob);
signal_prob_info.PixelDimensions = signal_prob_info.PixelDimensions(1:length(size(signal_prob))); % fixed in case only one signal IC is selected

potential_signal_prob = all_prob(:,:,:,signal_idx);
potential_signal_prob_info = noise_prob_info;
potential_signal_prob_info.ImageSize = size(potential_signal_prob);

niftiwrite(potential_signal_prob, 'PotentialSignalICs', potential_signal_prob_info, 'Compressed', true) % just the ICs that were considered
niftiwrite(noise_prob, 'NoiseICs', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob, 'SignalICs', signal_prob_info, 'Compressed', true)

% Now compress noise and signal, with only max probabilities per voxel
[~, ~] = max(all_prob, [], 4); 
[noise_prob_1D, noise_prob_ind] = max(all_prob(:,:,:,noise_ICs), [], 4); 
[signal_prob_1D, signal_prob_ind] = max(all_prob(:,:,:,signal_ICs), [], 4); 

% Update headers for the 3D images
noise_prob_info.ImageSize = noise_prob_info.ImageSize(1:end-1);
noise_prob_info.PixelDimensions = noise_prob_info.PixelDimensions(1:end-1);
signal_prob_info = noise_prob_info;

% Calculate out a signal to noise ratio mask overlap - but this might not
% be very helpful
signal_noise_ratio_IC_overlap = signal_prob_1D ./ (noise_prob_1D+signal_prob_1D);
signal_noise_ratio_IC_overlap(isnan(signal_noise_ratio_IC_overlap)) = 0.5; % convert if/when there is no noise or signal high prob value to 0.5 (will become 0)
signal_noise_ratio_IC_overlap = signal_noise_ratio_IC_overlap - 0.5; % shift 0.5 to 0, to center it so <0 is more nois & >0 is more signal

% perhaps a better method is threshold at 0.95 both signal and noise, make
% noise be noise - signal (to get noise sections alone), set to -1 for
% noise only sections, then combine signal and noise masks together.
thresholded_noise = noise_prob_1D > 0.949;
thresholded_signal = signal_prob_1D > 0.949;
noise_alone = thresholded_noise .* ~thresholded_signal;
signal_and_noise_overlap = thresholded_signal - noise_alone; % should just be 1 for signal, and -1 for noise

% write out the masks of highest probability
niftiwrite(noise_prob_1D, 'NoiseICOverlap', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob_1D, 'SignalICOverlap', signal_prob_info, 'Compressed', true)
%niftiwrite(signal_noise_ratio_IC_overlap, 'SignaltoNoiseICOverlap', signal_prob_info, 'Compressed', true) % <0 is more noise, >0 is more signal, 0 is either equivalent, or not high probability either way
%niftiwrite(cast(signal_noise_ratio_IC_overlap .* gm_bin, 'single'), 'SignaltoNoiseICOverlap_GM', signal_prob_info, 'Compressed', true) % <0 is more noise, >0 is more signal
niftiwrite(cast(signal_and_noise_overlap, 'single'), 'SignalandNoiseICOverlap', signal_prob_info, 'Compressed', true) % This is likely the most helpful one

% Write out the max noise ICs where we are in gray matter and noise is
% currently more represented, and vice versa. Can help in identifying
% misrepresnted signal and/or noise
int_prob_info = noise_prob_info;
int_prob_info.Datatype = 'int8';
niftiwrite(cast(noise_ICs(noise_prob_ind) .* gm_bin .* (signal_noise_ratio_IC_overlap < 0), 'int8'), 'Max_NoiseIC_GM', int_prob_info, 'Compressed', true)
niftiwrite(cast(signal_ICs(signal_prob_ind).*notgm_bin.*(signal_noise_ratio_IC_overlap>0), 'int8'), 'Max_SignalIC_NotGM', int_prob_info, 'Compressed', true)

% Export variables later use in next steps (e.g., fsl_regfilt)
% Notice there is "Manual" vs "Auto" prefix here, because these are the
% final manual picked ones
writematrix(Results.noise_ICs', 'manual_noise_dist_ICs.csv')
writematrix(Results.signal_ICs', 'manual_signal_dist_ICs.csv')
writematrix(Results.noise_indices, 'manual_noise_dist_indices.csv')
writematrix(Results.signal_indices, 'manual_signal_dist_indices.csv')

% and also the IC Checker table, normalized feature table, and compare
% cleaning estimates
writetable(Results.IC_checker_table, 'IC_manual_checker.csv')
writetable(Results_Auto.feature_relative_table, 'feature_manual_vals.csv')
writetable(Results_Auto.feature_table_norm, 'feature_manual_norms.csv')
writetable(Results.compare_cleaning, 'compare_manual_cleaning.csv','WriteRowNames', true)
writetable(Results.eval_table, 'eval_manual_table.csv')

% and then make noise components if you want to do aggressive denoising like in CONN (regression)
mixing_matrix = readmatrix([mel_fol, '/melodic_mix']);
manual_noise_dist_covariates = mixing_matrix(:, Results.noise_ICs');
writematrix(manual_noise_dist_covariates, 'manual_noise_dist_covariates.csv')
writematrix(manual_noise_dist_covariates, 'manual_noise_dist_covariates.txt', 'Delimiter', ' ')

% save a matrix of relevant variables if you want to examine later!
save('DecisionVariables_Manual.mat', 'output_dir', 'subject_id', 'session_id', ...
    'task_id', 'ICs', 'Data', 'Tables', 'Results_Auto', 'Results')

% And then move relevant files
movefile manual* ic_manual_selection
movefile *SignalIC* ic_manual_selection
movefile *NoiseIC* ic_manual_selection
movefile *IC_manual_checker.* ic_manual_selection
movefile *_manual_* ic_manual_selection
movefile *_Manual.* ic_manual_selection

if ~isfile('./ic_manual_selection/IC_manual_checker.csv')
    movefile(IC_manual_checker, "ic_manual_selection")
end

% OK, now run the manual CICADA denoising
prefix = [subject_id, '_', session_id, '_task-', task_id];
CICADA_tag = 'CICADA_manual_nonagg';
suffix = 'bold.nii.gz';
fsl_regfilt_command = ['fsl_regfilt -i ', task_dir, '/funcfile.nii.gz -f ', ...
    '"', '$(cat ', task_dir, '/ic_manual_selection/manual_noise_dist_ICs.csv)', '"', ...
    ' -d ', mel_fol, '/melodic_mix -m ', ...
    task_dir, '/funcmask.nii.gz -o ', prefix, '_', CICADA_tag, '_', suffix];
fprintf(['Running: ', fsl_regfilt_command, '\n'])
[~, ~] = call_fsl(fsl_regfilt_command);

% Now put things in cleaned_dir, and also get cleaned_filename
cleaned_dir = [task_dir, '/cleaned'];
movefile('sub*ses*task*CICADA*.nii.gz', cleaned_dir)
movefile('sub*ses*task**8p*.nii.gz', cleaned_dir)
movefile('sub*ses*task**9p*.nii.gz', cleaned_dir)

% get cleaned_filename
cleaned_file_info = dir([cleaned_dir, '/*CICADA*manual*.nii.gz']);
cleaned_file = [cleaned_file_info.folder, '/', cleaned_file_info.name];

end
