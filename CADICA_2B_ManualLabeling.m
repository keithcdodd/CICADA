function CADICA_2B_ManualLabeling(cadicafol, subject_id, session_id, taskname, IC_manual_checker)
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
% 
%
% IC_manual_checker (optional): IC_manual_checker.csv file - if not given, will look for it in ic_auto_selection
% folder

cd(cadicafol)
currsubjfol = ['sub-', subject_id];
cd(currsubjfol)
currsessfol = ['ses-', session_id];
cd(currsessfol)
cd(taskname)

if ~exist('IC_manual_checker', 'var')
    % If not given assume it is in its normal spot and look for it
    if isfile('./ic_auto_selection/IC_manual_checker.csv')
        IC_manual_checker = './ic_auto_selection/IC_manual_checker.csv';
    elseif isfile('./ic_manual_selection/IC_manual_checker.csv')
            IC_manual_checker = './ic_manual_selection/IC_manual_checker.csv';
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
if ~isfile('./ic_auto_selection/DecisionVariables_Auto.mat')
    fprintf('Cannot find ic_auto_selection/DecisionVariables_Auto.mat ... Exiting ...\n')
    return
end

load('./ic_auto_selection/DecisionVariables_Auto.mat')
auto_table = IC_checker_table;
IC_checker_table = manual_table; % update after loading it as desired

% remove old ic_manual_selection folder if it exists
% if isfolder('./ic_manual_selection')
%     rmdir './ic_manual_selection' 's'
% end

% Note what the changes are
IC_selection_changes = manual_table(manual_table.SignalLabel ~= auto_table.SignalLabel, :); % this will list the signal labels that manual has, that auto selection had the opposite on

writetable(IC_selection_changes, 'changed_signal_label_manual_ICs.csv')
ICs_signal_chance_order_table = array2table(signal_idx, 'VariableNames', {'ICs'});


% Now do the things
cut_off_table = IC_checker_table; 
if ismember('Classifications',cut_off_table.Properties.VariableNames)
    cut_off_table.Classifications = [];
end
manual_array = table2array(cut_off_table);
noise_indices(manual_array(manual_array(:,2) == 0,1)) = 1;
noise_indices(manual_array(manual_array(:,2) == 1,1)) = 0;
noise_indices = logical(noise_indices);

% Create final noise and signal ICs and indices
noise_ICs = ICs(noise_indices);
signal_indices = logical(~noise_indices);
signal_ICs = ICs(signal_indices);


% Make it easy to compare groups
close_ICs = int32(manual_array(:,1));
eval_names = ['IC Number', 'GM', 'Edge', 'Subependymal', 'CSF', 'Suscept', 'Outbrain', 'low freq' 'BOLDfreq', 'HighFreq', 'Spike_Corr', 'Move _Corr', 'Clustering_Prop', condition_names];
eval_table = array2table([ICs', feature_array], 'VariableNames', eval_names);
eval_table_close_ICs = eval_table(close_ICs, :);

% Compare before and after selection
before = sum(feature_array .* IC_exp_var, 1);
IC_exp_var_signal = IC_exp_var(signal_ICs);
IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
after = sum(feature_array(signal_ICs,:) .* IC_exp_var_signal, 1);
compare_cleaning = array2table([before', after'], 'RowNames', eval_names(2:end), 'VariableNames', {'Before', 'After'});

% write out compare_cleaning table, and the feature_norm table again
writetable(feature_table_norm, 'feature_manual_norms.csv') % this will be the same as the auto, but convenient to have in here too
writetable(compare_cleaning, 'compare_manual_cleaning.csv')
writetable(ICs_signal_chance_order_table, 'ICs_in_signal_chance_order.csv')

% Estimate effective degrees of freedom (somewhat) from CADICA alone via proportions:
percent_variance_kept = sum(IC_exp_var(signal_ICs)); % probably want > 0.1 (10%)
num_ICs_kept = length(signal_ICs);
num_ICs_total = length(ICs);

% Create images to review signal and noise components from
% probabilities ICs
all_prob = niftiread('./melodic/ICprobabilities.nii.gz');
all_clusterthresh = niftiread('./melodic/ICclusterthresh_zstat.nii.gz');

% higher probability clustered & Highest prob only clustered
all_prob_clustered = cast(all_prob .* (all_prob > 0.67) .* (abs(all_clusterthresh) > 3), 'single'); % includes probability values
all_highest_prob_clustered = cast(all_prob .* (all_prob > 0.95) .* (abs(all_clusterthresh) > 3), 'single'); % higher cut off

gm_prob = niftiread('./region_masks/GM_prop_final.nii.gz');
notgm_prob = niftiread('./region_masks/NotGM_prop_final.nii.gz');
gm_bin = gm_prob > 0.75;
notgm_bin = notgm_prob > 0.75;
noise_prob = all_prob_clustered(:,:,:,noise_ICs);
noise_prob_info = niftiinfo('./melodic/ICprobabilities.nii.gz');
signal_prob = all_prob_clustered(:,:,:,signal_ICs);

% grab noise probabilities, update header, and write, same for
% signal
noise_prob_info.ImageSize = size(noise_prob);
noise_prob_info.DisplayIntensityRange = [0.9 1];
noise_prob_info.Datatype = 'single';

signal_prob_info = noise_prob_info;
signal_prob_info.ImageSize = size(signal_prob);

potential_signal_prob = all_prob_clustered(:,:,:,potential_signal_ICs);
potential_signal_prob_info = noise_prob_info;
potential_signal_prob_info.ImageSize = size(potential_signal_prob);

niftiwrite(potential_signal_prob, 'PotentialSignalICs', potential_signal_prob_info, 'Compressed', true)
niftiwrite(noise_prob, 'NoiseICs', noise_prob_info, 'Compressed', true)
% I need to code something that can handle if it is only one signal, and
% therefore only 3D instead of expected 4D
niftiwrite(signal_prob, 'SignalICs', signal_prob_info, 'Compressed', true)

% Now compress noise and signal, with only highest probabilities
[all_prob_1D, all_prob_ind] = max(all_prob, [], 4); 
[noise_prob_1D, noise_prob_ind] = max(all_prob(:,:,:,noise_ICs), [], 4); 
[signal_prob_1D, signal_prob_ind] = max(all_prob(:,:,:,signal_ICs), [], 4); 

% Update headers for the 3D images
noise_prob_info.ImageSize = noise_prob_info.ImageSize(1:end-1);
noise_prob_info.PixelDimensions = noise_prob_info.PixelDimensions(1:end-1);
signal_prob_info = noise_prob_info;

% Calculate out a signal to noise ratio mask overlap
signal_noise_ratio_IC_overlap = signal_prob_1D ./ (noise_prob_1D+signal_prob_1D);
signal_noise_ratio_IC_overlap(isnan(signal_noise_ratio_IC_overlap)) = 0.5; % convert if/when there is no noise or signal high prob value to 0.5 (will become 0)
signal_noise_ratio_IC_overlap = signal_noise_ratio_IC_overlap - 0.5; % shift 0.5 to 0, to center it so <0 is more noise & >0 is more signal

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
niftiwrite(signal_noise_ratio_IC_overlap, 'SignaltoNoiseICOverlap', signal_prob_info, 'Compressed', true) % <0 is more noise, >0 is more signal, 0 is either equivalent, or not high probability either way
niftiwrite(cast(signal_noise_ratio_IC_overlap .* gm_bin, 'single'), 'SignaltoNoiseICOverlap_GM', signal_prob_info, 'Compressed', true) % <0 is more noise, >0 is more signal
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
writematrix(noise_ICs, 'manual_noise_dist_ICs.csv')
writematrix(signal_ICs, 'manual_signal_dist_ICs.csv')
writematrix(noise_indices, 'manual_noise_dist_indices.csv')
writematrix(signal_indices, 'manual_signal_dist_indices.csv')

% and then make noise components if you want to do aggressive denoising like in CONN (regression)
mixing_matrix = dlmread('./melodic/melodic_mix');
manual_noise_dist_covariates = mixing_matrix(:, noise_ICs);
writematrix(manual_noise_dist_covariates, 'manual_noise_dist_covariates.csv')
writematrix(manual_noise_dist_covariates, 'manual_noise_dist_covariates.txt', 'Delimiter', ' ')

% save a matrix of relevant variables if you want to examine later!
save('DecisionVariables_Manual.mat', 'cadicafol', 'subject_id', 'session_id', 'taskname', 'ICs', 'classification_table',...
    'eval_table', 'signal_ICs', 'noise_ICs', 'signal_checker_table', 'IC_checker_table', 'repetitiontime', 'numvolumes',...
    'feature_array', 'feature_table', 'IC_exp_var', 'compare_cleaning', 'close_ICs', 'noise_indices', 'HighSig_ICs', 'signal_idx',...
    'Highnoise_ICs', 'averages_table', 'std_table', 'classification_table', 'percent_variance_kept', ...
    'IC_exp_var', 'num_ICs_kept', 'num_ICs_total', 'Highconditioncorr_ICs', 'feature_table_norm')

movefile manual* ic_manual_selection
movefile *SignalIC* ic_manual_selection
movefile *NoiseIC* ic_manual_selection
movefile *Manual* ic_manual_selection
movefile *manual_ICs.* ic_manual_selection
if ~isfile('./ic_manual_selection/IC_manual_checker.csv')
    movefile(IC_manual_checker, "ic_manual_selection")
end

end
