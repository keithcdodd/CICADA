function T_results = rank_ICs_by_group_NSP(home_dir, IC_mel_data, group_funcmask_file, subj_func_files, subj_notgm_files, TR, FD_cell, DVARS_cell, hrf_general_power, hrf_task_power, hrf_conditions)
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
max_condition_corr_idx = zeros(nSubjects, nICs);
power_spectra_all = {};  % Cell array: {subject}(T/2+1 × nICs)
fd_corr_subjectwise = NaN(nSubjects, nICs);
dvars_corr_subjectwise = NaN(nSubjects, nICs);
if ~isempty(hrf_task_power)
    condition_corrs = zeros(nICs, size(hrf_task_power, 2), nSubjects); % ICs x conditions x nSubjects
else
    condition_corrs = [];
end

%% Loop over subjects
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
    power_overlap_general = sum(hrf_general_power .* ps', 2); % 1 x nICs

end

% Finish with hrf_task_power overlap
% idea here is we take the overlap of the condition corr that best matches
% the given IC across all subjects
if ~isempty(condition_corrs)
    % need to figure out the mode of index of condition_corrs per IC
    max_condition_corrs_indices = zeros(size(nICs, nSubjects));
    mode_max_condition_corrs_indices = zeros(size(nICs));
    max_conditioncorr_powerspectra = zeros(size(hrf_task_power, 1), nICs);
    for i = 1:nICs
        for s = 1:nSubjects
            [~,curr_I] = max(condition_corrs(i, :, s), [], 2);
            max_condition_corrs_indices(i, s) = curr_I; % which condition was most correlated to the IC timeseries for each subject
        end
        mode_max_condition_corrs_indices(i) = mode(max_condition_corrs_indices(i, :)); % mode is the most represented condition across all subjects for each IC
        max_conditioncorr_powerspectra(:,i) = hrf_task_power(:, mode_max_condition_corrs_indices(i))'; % powerspectra to use per IC for overlap comparison!
    end

    power_overlap_all_hrfs = sum(max_conditioncorr_powerspectra' .* ps', 2); % I want this to output the subject average task overlap per IC (1 x nIC)

    % OK, now take max value between this and power_overlap_general
    [power_overlap_max, power_overlap_max_idx] = max([power_overlap_all_hrfs(:), power_overlap_general(:)], [], 2); % 1 x nICs, idx is 1 for conditions, 2 for general
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
'VariableNames', {'IC', 'Signal Label Guess', 'Group_NSP_norm', 'GM_Overlap_Norm', 'Power_Overlap_Norm', 'Smoothness_Norm', 'Power_Overlap_Condition', 'FD_Correlation_Norm', 'DVARS_Correlation_Norm'});
T_results = T_results(sorted_idx, :);


% %% Plot average power spectra
% freqs = fs * (0:floor(T/2)) / T;
% figure;
% hold on;
% colors = parula(nICs);
% for i = 1:nICs
%     plot(freqs, mean_power_spectra(:, sorted_idx(i)), ...
%         'Color', colors(sorted_idx(i), :), 'LineWidth', 1.2);
% end
% xlabel('Frequency (Hz)');
% ylabel('Power');
% title('Average Power Spectrum per IC');
% legend(arrayfun(@(i) sprintf('IC %d', i), sorted_idx(1:nICs), 'UniformOutput', false), ...
%     'Location', 'northeastoutside');
% grid on;
% hold off;

%% Optional: Save per-IC summary plots (uncomment to use)
group_ica_eval_dir = [home_dir, '/ica_evals'];
if ~isfolder(group_ica_eval_dir)
    mkdir(group_ica_eval_dir)
end
for i = 1:nICs
    ic_idx = sorted_idx(i);
    figure('Visible','off');
    subplot(2,1,1);
    plot(freqs, mean_power_spectra(:, ic_idx), 'b', 'LineWidth', 1.5);
    xlabel('Frequency (Hz)');
    ylabel('Power');
    title(sprintf('IC %d - Power Spectrum', ic_idx));

    subplot(2,1,2);
    bar([1 2 3], [mean_gm_overlap_norm(ic_idx), mean_abs_fd_corr_norm(ic_idx), mean_abs_dvars_corr_norm(ic_idx)]);
    text(1:3,[mean_gm_overlap_norm(ic_idx), mean_abs_fd_corr_norm(ic_idx), mean_abs_dvars_corr_norm(ic_idx)],num2str([round(mean_gm_overlap_norm(ic_idx), 2), round(mean_abs_fd_corr_norm(ic_idx), 2), round(mean_abs_dvars_corr_norm(ic_idx), 2)]'),'vert','bottom','horiz','center'); 
    box off
    set(gca, 'XTickLabel', {'GM Overlap Scaled', 'FD Corr', 'DVARS Corr'});
    %ylim([-1.1 1.1]);
    title(sprintf('IC %d - Overlap and Motion Corrs', ic_idx));

    saveas(gcf, [group_ica_eval_dir, '/GM_', num2str(i), '_IC_', num2str(ic_idx), '_eval.png']);
    close;
end

writetable(T_results, [group_ica_eval_dir, '/IC_eval_results.csv'], 'Delimiter',',')
disp(T_results)
end
