%% Updates to signal labeling by manual checking 
% This script is to allow you to examine the auto labels, and quickly
% manually adjust them as you see fit. Follow directions below. We suggest
% using the Griffanti 2017 Hand classification of fMRI ICA noise components
% paper as a reference. This needs to be run task by task and therefore is
% the longest step. The point is that this whole process gives the
% same/similar effect of manual ICA denoising, but it makes it much faster
% and easier while providing more information.


%% 1. Run the Small Code Block (after making sure variables are correct)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars
%%%%% Variables to change/check
cadicafol='/Volumes/VectoTec_VectoTech_Media_Rapid/DMXBA/CADICA_Updated';
% similar syntax as in the 1_CADICA_MasksandICAs
subject = '107';
session = '03';
task = 'rest';
task_names = {}; %empty if rest

cd(cadicafol)
cd(['sub-', subject])
cd(['ses-', session])
cd(task)

load('DecisionVariables_Auto.mat') % load the relevant materials
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 2. Open IC_checker_updated_table variable, and update signal labels as you see fit. 
% Suggestion: Use table and classification table, alongside the MELODIC
% report and the melodic_IC.nii.gz file overlayed on anatfile in FSLEYES to
% help make decisions

%% 3. After updating IC_checker_updated_table, run the following code block below to finish up!

close_IC_checker_updated_array = table2array(IC_checker_updated_table);
noise_indices(close_IC_checker_updated_array(close_IC_checker_updated_array(:,2) == 0,1)) = 1;
noise_indices(close_IC_checker_updated_array(close_IC_checker_updated_array(:,2) == 1,1)) = 0;
noise_indices = logical(noise_indices);

% Create final noise and signal ICs and indices, if GM_prop>0.5,
% label as signal regardless
% check close ICs and update noise indices
noise_ICs = ICs(noise_indices);
signal_indices = logical(~noise_indices);
signal_ICs = ICs(signal_indices);


% Make it easy to compare groups
close_ICs = int32(close_IC_checker_updated_array(:,1));
eval_names = ['IC Number', 'GM', 'Edge', 'Subependymal', 'CSF', 'Suscept', 'Outbrain', 'low freq' 'BOLDfreq', 'HighFreq', 'Spike_Corr', 'Move _Corr', 'Clustering_Prop', task_names];
eval_table = array2table([ICs', feature_array], 'VariableNames', eval_names);
eval_table_close_ICs = eval_table(close_ICs, :);

% Compare before and after selection
before = sum(feature_array .* IC_exp_var);
IC_exp_var_signal = IC_exp_var(signal_ICs);
IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
after = sum(feature_array(signal_ICs,:) .* IC_exp_var_signal);
compare_cleaning = array2table([before', after'], 'RowNames', eval_names(2:end), 'VariableNames', {'Before', 'After'});

% Estimate effective degrees of freedom (somewhat) from CADICA alone via proportions:
percent_variance_kept = sum(IC_exp_var(signal_ICs)); % probably want > 0.1 (10%)

% Create images to review signal and noise components from
% probabilities ICs
noise_prob = niftiread('./melodic/ICprobabilities.nii.gz');
noise_prob_info = niftiinfo('./melodic/ICprobabilities.nii.gz');
signal_prob = noise_prob; signal_prob_info = noise_prob_info;

% Grab noise probabilities, update header, and write, same for signal
noise_prob = noise_prob(:,:,:,noise_ICs);
noise_prob_info.ImageSize = size(noise_prob);
noise_prob_info.DisplayIntensityRange = [0.9 1];
signal_prob = signal_prob(:,:,:, signal_ICs);
signal_prob_info = noise_prob_info;
signal_prob_info.ImageSize = size(signal_prob);

niftiwrite(noise_prob, 'NoiseICs', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob, 'SignalICs', signal_prob_info, 'Compressed', true)

% Now binarize and compress noise and signal
noise_prob(noise_prob>0.899) = 1; noise_prob(noise_prob<0.9) = 0;
noise_prob = sum(noise_prob, 4); noise_prob(noise_prob>0.99) = 1;

% Grab, binarize and sum signal_probs
signal_prob(signal_prob>0.899) = 1; signal_prob(signal_prob<0.9) = 0;
signal_prob = sum(signal_prob, 4); signal_prob(signal_prob>0.99) = 1;

% Update headers for the 3D images
noise_prob_info.ImageSize = noise_prob_info.ImageSize(1:end-1);
noise_prob_info.PixelDimensions = noise_prob_info.PixelDimensions(1:end-1);
signal_prob_info = noise_prob_info;

% write out the masks
niftiwrite(noise_prob, 'NoiseICOverlap', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob, 'SignalICOverlap', signal_prob_info, 'Compressed', true)

% (5c) Export variables later use in next steps (e.g., fsl_regfilt)
% Notice there is "Manual" vs "Auto" prefix here, because these are the
% final manual picked ones
writematrix(noise_ICs, 'Manual_Noise_dist_ICs.csv')
writematrix(signal_ICs, 'Manual_Signal_dist_ICs.csv')
writematrix(noise_indices, 'Manual_Noise_dist_indices.csv')
writematrix(signal_indices, 'Manual_Signal_dist_indices.csv')

% and then make noise components if you want to do aggressive denoising like in CONN (regression)
mixing_matrix = dlmread('./melodic/melodic_mix');
Manual_noise_dist_covariates = mixing_matrix(:, noise_ICs);
save('Manual_CADICA_Noise_dist.mat', 'Manual_noise_dist_covariates')

% save a matrix of relevant variables if you want to examine later!
save('DecisionVariables_Manual.mat', 'ICs', 'classification_table',...
    'eval_table', 'signal_ICs', 'noise_ICs', 'signal_checker_table', 'IC_checker_updated_table', ...
    'feature_array', 'feature_table', 'IC_exp_var', 'compare_cleaning', 'close_ICs', 'noise_indices', 'HighSig_ICs', ...
    'Highnoise_ICs', 'averages_table', 'std_table', 'classification_table', 'percent_variance_kept', 'Hightaskcorr_ICs', 'feature_table_norm')


