function [confusion_matrix_table, predictive_matrix_table, noise_selection_table, signal_selection_table] = confusion_matrix(task_dir)
% Generate confusion matrix comparing manual to other methods
% task_dir is the CICADA task directory for a subject

% First, look for ic_manual_selection folder with manual_noise_dist_indices.csv and
% manual_signal_dist_indices.csv

if ~isfile([task_dir, '/ic_manual_selection/manual_noise_dist_indices.csv']) || ~isfile([task_dir, '/ic_manual_selection/manual_signal_dist_indices.csv'])
    fprintf('Cannot find both the noise and signal indices from manual CICADA.\n')
    return
end

manual_noise_indices = readmatrix([task_dir, '/ic_manual_selection/manual_noise_dist_indices.csv']);
manual_signal_indices = readmatrix([task_dir, '/ic_manual_selection/manual_signal_dist_indices.csv']);

ICs = 1:length(manual_noise_indices);
base_indices = zeros(length(ICs),1);

noise_selection_labels = "Manual";
signal_selection_labels = "Manual";
noise_selection_array = manual_noise_indices;
signal_selection_array = manual_signal_indices;

% OK, now step through looking for each other folder of potential interest
% to create comparison table
% First, AutoCICADA:
curr_noise_ICs = [task_dir, '/ic_auto_selection/auto_noise_dist_ICs.csv'];
curr_label = "AutoCICADA";
if isfile(curr_noise_ICs)
    % Now we can convert noise labels to indices
    curr_noise_ICs = readmatrix(curr_noise_ICs);
    curr_noise_indices = base_indices; 
    curr_noise_indices(curr_noise_ICs) = 1;
    curr_signal_indices = ~curr_noise_indices;

    % update everything
    noise_selection_labels = [noise_selection_labels, curr_label];
    signal_selection_labels = [signal_selection_labels, curr_label];

    noise_selection_array = [noise_selection_array, curr_noise_indices];
    signal_selection_array = [signal_selection_array, curr_signal_indices];
end

% FIX
curr_noise_ICs = [task_dir, '/fix/fix_noise_ICs.csv'];
curr_label = "Fix";
if isfile(curr_noise_ICs)
    % Now we can convert noise labels to indices
    curr_noise_ICs = readmatrix(curr_noise_ICs);
    curr_noise_indices = base_indices; 
    curr_noise_indices(curr_noise_ICs) = 1;
    curr_signal_indices = ~curr_noise_indices;

    % update everything
    noise_selection_labels = [noise_selection_labels, curr_label];
    signal_selection_labels = [signal_selection_labels, curr_label];

    noise_selection_array = [noise_selection_array, curr_noise_indices];
    signal_selection_array = [signal_selection_array, curr_signal_indices];
end

% AROMA
curr_noise_ICs = [task_dir, '/ic_aroma_selection/classified_motion_ICs.txt'];
curr_label = "Aroma";
if isfile(curr_noise_ICs)
    % Now we can convert noise labels to indices
    curr_noise_ICs = readmatrix(curr_noise_ICs);
    curr_noise_indices = base_indices; 
    curr_noise_indices(curr_noise_ICs) = 1;
    curr_signal_indices = ~curr_noise_indices;

    % update everything
    noise_selection_labels = [noise_selection_labels, curr_label];
    signal_selection_labels = [signal_selection_labels, curr_label];

    noise_selection_array = [noise_selection_array, curr_noise_indices];
    signal_selection_array = [signal_selection_array, curr_signal_indices];
end

% ALT
curr_noise_ICs = [task_dir, '/ic_alt_selection/noise_labels.txt'];
curr_label = "Alt";
if isfile(curr_noise_ICs)
    % Now we can convert noise labels to indices
    curr_noise_ICs = readmatrix(curr_noise_ICs);
    curr_noise_indices = base_indices; 
    curr_noise_indices(curr_noise_ICs) = 1;
    curr_signal_indices = ~curr_noise_indices;

    % update everything
    noise_selection_labels = [noise_selection_labels, curr_label];
    signal_selection_labels = [signal_selection_labels, curr_label];

    noise_selection_array = [noise_selection_array, curr_noise_indices];
    signal_selection_array = [signal_selection_array, curr_signal_indices];
end


% Create confusion matrices:
true_noise_array = zeros(size(1, width(noise_selection_array)));
true_signal_array = true_noise_array;
false_noise_array = true_noise_array;
false_signal_array = true_noise_array;
for idx = 1:width(noise_selection_array)
    % first slot should always be 100% because it is comparing manual to
    % manual
    true_noise_array(idx) = sum(noise_selection_array(:,1) == 1 & noise_selection_array(:,idx) == 1) ./ length(ICs);
    true_signal_array(idx) = sum(signal_selection_array(:,1) == 1 & signal_selection_array(:,idx) == 1) ./ length(ICs);

    false_noise_array(idx) = sum(noise_selection_array(:,1) == 0 & noise_selection_array(:,idx) == 1) ./ length(ICs);
    false_signal_array(idx) = sum(signal_selection_array(:,1) == 0 & signal_selection_array(:,idx) == 1) ./ length(ICs);
end
confusion_matrix_extended = [true_noise_array; false_noise_array; true_signal_array; false_signal_array]; % lets you compare to manual percentages as first column
confusion_matrix = [true_noise_array(2:end); false_noise_array(2:end); true_signal_array(2:end); false_signal_array(2:end)];


sensitivity_matrix = confusion_matrix(1,:) ./ (confusion_matrix(1,:) + confusion_matrix(4,:)); % if it is noise, do we capture it?
specificity_matrix = confusion_matrix(3,:) ./ (confusion_matrix(3,:) + confusion_matrix(2,:)); % if it is signal, do we capture it?
precision_matrix = confusion_matrix(1,:) ./ (confusion_matrix(1,:) + confusion_matrix(2,:)); % same as positive predictive value. If we label it as noise, how sure are we that is correct?
neg_pred_value_matrix = confusion_matrix(3,:) ./ (confusion_matrix(3,:) + confusion_matrix(4,:)); % negative predictive value. If we label it as signal, how sure are we that is correct?
accuracy_matrix = ( confusion_matrix(1,:) + confusion_matrix(3,:) ) ./ sum(confusion_matrix); % overall, how accurate are we with noise vs signal labels?

predictive_matrix = [sensitivity_matrix; precision_matrix; specificity_matrix; neg_pred_value_matrix; accuracy_matrix];

% now put everything together into tables!
noise_selection_table = array2table(noise_selection_array, 'VariableNames', noise_selection_labels);
signal_selection_table = array2table(signal_selection_array, 'VariableNames', signal_selection_labels);
% Note, positive means capturing noise IC, negative means capturing signal
% IC. This follows convention of AROMA paper.
confusion_matrix_table = array2table(confusion_matrix_extended, 'VariableNames', noise_selection_labels, 'RowNames', ["True_Noise", "False_Noise", "True_Signal", "False_Signal"]);
predictive_matrix_table = array2table(predictive_matrix, 'VariableNames', noise_selection_labels(2:end), 'RowNames', ["Noise Sensitivity", "Noise_Predictive_Value", "Signal Specificity", "Signal_Predictive_Value", "Accuracy"]);

writetable(noise_selection_table, 'noise_selection_table.csv')
writetable(signal_selection_table, 'signal_selection_table.csv')
writetable(confusion_matrix_table, 'confusion_matrix_table.csv', 'WriteRowNames',true)
writetable(predictive_matrix_table, 'predictive_matrix_table.csv', 'WriteRowNames',true)


end