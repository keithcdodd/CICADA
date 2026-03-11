function T_results = rank_ICs_by_group_NSP(output_dir, IC_mel_data, group_funcmask_file, subj_func_files, subj_notgm_files, TR, FD_cell, DVARS_cell, hrf_general_power, hrf_task_power, hrf_conditions)
% By ranking group ICs by GM overlap and computing correlation between
% the magnitude of the temporal derivative of each IC timecourse and FD/DVARS,
% this function helps distinguish signal from noise components for a user
% (see GM overlap, power spectra, and FD and DVARS correlation).
%
% Inputs:
% - IC_file: path to group-level melodic_IC.nii.gz
% - subj_func_files: cell array of subject fMRI .nii.gz files
% - subj_gm_files: cell array of subject GM probability maps
% - TR: repetition time in seconds
% - FD_cell: cell array of FD vectors (T-1 elements each)
% - DVARS_cell: cell array of DVARS vectors (T-1 elements each)
% - hrf_general_power: to sum overlap to powerspectra (resting & task)
% - hrf_task_power: to sum overlap to powerspectra (task only)
% - hrf_conditions (table): blocks convolved with hrf for each condition 
%
% Output:
% - T_results: table of ICs sorted by GM overlap with FD and DVARS correlation

nSubjects = length(subj_func_files);
if ~isempty(hrf_conditions)
    hrf_condition_names = hrf_conditions.Properties.VariableNames;
    hrf_conditions = table2array(hrf_conditions); % convert to array
end

%% Load Group IC Map
ICs = IC_mel_data;  % X × Y × Z × N_ICs
[nx, ny, nz, nICs] = size(ICs);
IC_2D = reshape(ICs, [], nICs);   % Voxel × ICs
IC_2D = IC_2D ./ max(abs(IC_2D), [], 1);  % normalize each IC

%% load group functional mask
group_funcmask = niftiread(group_funcmask_file);
group_funcmask_info = niftiinfo(group_funcmask_file);
voxel_size = group_funcmask_info.PixelDimensions;

%% Initialize metrics
gm_overlap_all = zeros(nSubjects, nICs);
fwhm_avg_all = zeros(nSubjects, nICs);
power_spectra_all = {};  % Cell array: {subject}(T/2+1 × nICs)
fd_corr_subjectwise = NaN(nSubjects, nICs);
dvars_corr_subjectwise = NaN(nSubjects, nICs);
if ~isempty(hrf_task_power)
    condition_corrs = zeros(nICs, size(hrf_task_power, 2), nSubjects); % ICs x conditions x nSubjects
else
    condition_corrs = [];
end

%% Loop over subjects
power_overlap_general_all = zeros(nSubjects, nICs);
for s = 1:nSubjects
    fprintf('Processing subject %d/%d\n', s, nSubjects);

    % Load functional data
    func_data = niftiread(subj_func_files{s}); % X × Y × Z × T
    T = size(func_data, 4);
    func_2D = reshape(func_data, [], T);  % Voxel × Time

    % Mask: valid voxels
    valid_vox = all(~isnan(func_2D), 2) & any(func_2D ~= 0, 2);
    func_2D_valid = func_2D(valid_vox, :);
    IC_valid = IC_2D(valid_vox, :);

    %% Dual Regression Stage 1
    timecourses = (pinv(IC_valid) * func_2D_valid)';  % T × N_ICs

    %% Dual Regression Stage 2
    recon_data = func_2D_valid * pinv(timecourses)';  % Voxel × N_ICs
    subj_IC_maps = zeros(size(func_2D, 1), nICs);
    subj_IC_maps(valid_vox, :) = recon_data;

    

    %% Load GM map
    notGM = niftiread(subj_notgm_files{s});
    GM = ~notGM .* group_funcmask;
    GM_1D = reshape(GM, [], 1);

    %% GM overlap & smoothness
    subj_IC_maps = subj_IC_maps ./ max(abs(subj_IC_maps), [], 1);  % normalize
    for i = 1:nICs
        gm_overlap_all(s, i) = sum(subj_IC_maps(:, i) .* GM_1D) / sum(GM_1D);
        % compute smoothness per IC with FSL smoothest
        fwhm_avg_all(s,i) = estimate_ic_smoothness(subj_IC_maps(:,i), group_funcmask_file, [nx, ny, nz], voxel_size);
    end

    %% Power spectrum + FD/DVARS correlations
    fs = 1 / TR;
    freqs = fs * (0:floor(T/2)) / T;
    ps = zeros(length(freqs), nICs);

    for i = 1:nICs
        ts = detrend(timecourses(:, i));
        Y = fft(ts);
        P2 = abs(Y / T).^2;
        P1 = P2(1:floor(T/2)+1);
        P1(2:end-1) = 2 * P1(2:end-1);
        P1_norm = P1 ./ trapz(P1);
        ps(:, i) = P1_norm;
        

        %% Motion correlation using timecourse derivative
        ts_deriv = abs(diff(ts));  % T-1 × 1

        fd = FD_cell{s};
        fd = fd(2:end); % remove first NAN timepoint
        dvars = DVARS_cell{s};
        dvars = dvars(2:end); % remove first NAN timepoint 

        if length(fd) ~= T-1 || length(dvars) ~= T-1
            error('FD and DVARS must have length T-1 for subject %d (T = %d)', s, T);
        end

        fd_corr_subjectwise(s, i) = corr(ts_deriv, fd(:), 'Rows', 'complete');
        dvars_corr_subjectwise(s, i) = corr(ts_deriv, dvars(:), 'Rows', 'complete');

        % hrf task (if exists) and general:
        % load HRF info and get general and task overlap (if it exists)
        if ~isempty(condition_corrs)
            % is task data, find overlap of most correlated 
        
            for m = 1:size(condition_corrs, 2)
                condition_corrs(i,m,s) = corr(hrf_conditions(:,m), ts).^2; % ICs x conditions x nSubjects
            end
        end
    end

    power_spectra_all{s} = ps;
    % Now do the easier hrf general that definitely exists
    power_overlap_general_all(s, :) = sum(hrf_general_power .* ps', 2)';

end

power_overlap_general = mean(power_overlap_general_all, 1); % I want this to output the subject average general hrf overlap per IC (1 x nIC)

% Finish with hrf_task_power overlap
% idea here is we take the overlap of the condition corr that best matches
% the given IC across all subjects
if ~isempty(condition_corrs)
    % need to figure out the mode of index of condition_corrs per IC
    max_condition_corrs_indices = zeros(nICs, nSubjects);
    mode_max_condition_corrs_indices = zeros(nICs, 1);
    max_conditioncorr_powerspectra = zeros(size(hrf_task_power, 1), nICs);
    for i = 1:nICs
        for s = 1:nSubjects
            [~,curr_I] = max(condition_corrs(i, :, s), [], 2);
            max_condition_corrs_indices(i, s) = curr_I; % which condition was most correlated to the IC timeseries for each subject
        end
        mode_max_condition_corrs_indices(i) = mode(max_condition_corrs_indices(i, :)); % mode is the most represented condition across all subjects for each IC
        max_conditioncorr_powerspectra(:,i) = hrf_task_power(:, mode_max_condition_corrs_indices(i))'; % powerspectra to use per IC for overlap comparison!
    end

    power_overlap_task_all = zeros(nSubjects, nICs);

    for s = 1:nSubjects
        curr_ps = power_spectra_all{s};  % freq x nICs
        power_overlap_task_all(s, :) = sum(max_conditioncorr_powerspectra .* curr_ps, 1);
    end
    
    power_overlap_all_hrfs = mean(power_overlap_task_all, 1);  % 1 x nICs, mean over subjects

    % OK, now take max value between this and power_overlap_general
    [power_overlap_max, power_overlap_max_idx] = max([power_overlap_general(:), power_overlap_all_hrfs(:)], [], 2); % 1 x nICs, idx is 1 for general, 2 for tasks
    power_overlap_max = power_overlap_max(:); % make sure it is a column
else
    power_overlap_all_hrfs = [];
    power_overlap_max = power_overlap_general(:); % make sure it is a column
    power_overlap_max_idx = ones(size(power_overlap_max)); % all just based on general HRF
end

power_overlap_condition = cell(size(power_overlap_max_idx));
for idx = 1:length(power_overlap_max_idx)
    if power_overlap_max_idx(idx) == 1
        power_overlap_condition{idx} = 'General';
    else
        % we are working with a task condition then
        power_overlap_condition{idx} = hrf_condition_names{mode_max_condition_corrs_indices(idx)};
    end
end


%% Average metrics across subjects
mean_gm_overlap = mean(gm_overlap_all, 1);
mean_power_overlap = power_overlap_max;
mean_smoothness = mean(fwhm_avg_all, 1);
mean_power_spectra = mean(cat(3, power_spectra_all{:}), 3);  % freq × ICs
mean_abs_fd_corr = mean(abs(fd_corr_subjectwise), 1, 'omitnan');
mean_abs_dvars_corr = mean(abs(dvars_corr_subjectwise), 1, 'omitnan');

%% Group NSP (roughly)
gNSP = mean_gm_overlap(:) .* mean_power_overlap(:) .* mean_smoothness(:);

%% Rank by group NSP (GM, power overlap, smoothness)
[~, sorted_idx] = sort(gNSP, 'descend');

gNSP_norm = normalize(gNSP, 'range');
mean_gm_overlap_norm = normalize(mean_gm_overlap, 'range'); % easier to interpret
mean_power_overlap_norm = normalize(mean_power_overlap, 'range');
mean_smoothness_norm = normalize(mean_smoothness, 'range');
mean_abs_fd_corr_norm = normalize(mean_abs_fd_corr, 'range');
mean_abs_dvars_corr_norm = normalize(mean_abs_dvars_corr, 'range');

% kmeans cluster signal and not
[idx, ~] = kmeans(gNSP_norm(:), 3, 'Start', [0; median(gNSP_norm(:)); 1]);
signal_label = double(idx ~= 1)';

%% Output table
T_results = table((1:nICs)', signal_label', gNSP_norm(:), mean_gm_overlap_norm', mean_power_overlap_norm(:), mean_smoothness_norm', power_overlap_condition, mean_abs_fd_corr_norm', mean_abs_dvars_corr_norm', ...
'VariableNames', {'IC', 'Signal_Label', 'Group_NSP_norm', 'GM_Overlap_Norm', 'Power_Overlap_Norm', 'Smoothness_Norm', 'Power_Overlap_Condition', 'FD_Correlation_Norm', 'DVARS_Correlation_Norm'});
T_results = T_results(sorted_idx, :);

group_ica_eval_dir = fullfile(output_dir, 'ica_evals');
if ~isfolder(group_ica_eval_dir)
    [ok, msg, msgID] = mkdir(group_ica_eval_dir);
    if ~ok
        error('Failed to create folder "%s"\nmsg: %s\nmsgID: %s', ...
            group_ica_eval_dir, msg, msgID);
    end
end

for i = 1:nICs
    ic_idx = sorted_idx(i);

    f = figure('Visible','off', 'Color','w');

    subplot(2,1,1);
    plot(freqs, mean_power_spectra(:, ic_idx), 'b', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(sprintf('IC %d - Power Spectrum', ic_idx));
    box off

    subplot(2,1,2);
    vals = [mean_gm_overlap_norm(ic_idx), ...
            mean_abs_fd_corr_norm(ic_idx), ...
            mean_abs_dvars_corr_norm(ic_idx)];

    bar([1 2 3], vals);
    text(1:3, vals, compose('%.2f', vals), ...
        'VerticalAlignment', 'bottom', ...
        'HorizontalAlignment', 'center');
    box off
    set(gca, 'XTick', [1 2 3], ...
             'XTickLabel', {'GM Overlap Scaled', 'FD Corr', 'DVARS Corr'});
    title(sprintf('IC %d - Overlap and Motion Corrs', ic_idx));

    out_file = fullfile(group_ica_eval_dir, ...
        sprintf('GM_%d_IC_%d_eval.png', i, ic_idx));

    save_figure_robust(f, out_file, 300);
    close(f);
end

writetable(T_results, fullfile(group_ica_eval_dir, 'IC_eval_results.csv'), 'Delimiter', ',')
disp(T_results)
end
