function plot_qc(denoised_Edge_Edge_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_Outbrain_corr, ...
    denoised_WMCSF_WMCSF_corr, denoised_CSF_CSF_corr, denoised_NotGM_NotGM_corr, denoised_GM_GM_corr, denoised_Suscept_Suscept_corr, ...
    compare_Edge_Edge_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, compare_Outbrain_Outbrain_corr, ...
    compare_WMCSF_WMCSF_corr, compare_CSF_CSF_corr, compare_NotGM_NotGM_corr, compare_GM_GM_corr, compare_Suscept_Suscept_corr, ...
    orig_Edge_Edge_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, orig_Outbrain_Outbrain_corr, ...
    orig_WMCSF_WMCSF_corr, orig_CSF_CSF_corr, orig_NotGM_NotGM_corr, orig_GM_GM_corr, orig_Suscept_Suscept_corr, ...
    denoised_GM_mean, compare_GM_mean, orig_GM_mean, title_string, qc_plot_dest, denoised_name)

% Plot QC results! Add this in (will need to add in compare and orig corrs
% somehow)

% check if inputs are cell arrays instead of doubles. If so, convert to
% doubles


% There needs to be a decision on whether we are including GM mean signal
% (good for within subject comparison, worthless in group qc comparison)
% You can see if GM mean is empty for denoised

fprintf('Creating Figures\n')
if isempty(denoised_GM_mean)
    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 24, 10.8]);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 12, 5.4]);
    t = tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');
else
    fig = figure('Visible', 'off', 'Units', 'inches', 'Position', [1, 1, 24, 14.4]);
    set(fig, 'PaperUnits', 'inches', 'PaperPosition', [0, 0, 24, 14.4]);
    t = tiledlayout(4,3, 'Padding', 'compact', 'TileSpacing', 'compact');
end

lgd_font_size = 9;
axes_font_size = 11;

set(fig, 'DefaultLegendFontSize', lgd_font_size);
set(fig, 'DefaultLegendBox', 'off');
set(fig, 'DefaultAxesFontSize', axes_font_size);


%if isempty(denoised_GM_mean)
%	fig = figure('Position', [50, 200, 1700, 750], 'Visible', 'off');
%    t = tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');
%else
%    fig = figure('Position', [50, 50, 1700, 1000], 'Visible', 'off');
%    t = tiledlayout(4,3, 'Padding', 'compact', 'TileSpacing', 'compact');
%end

title(t, title_string, 'Interpreter', 'none')

% GM_Edge corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
colorord = get(gca, 'colororder'); % this has to be after the nexttile call because of some weird Matlab bug
hold on
title('Edge Correlation', 'Interpreter', 'none')
if ~isempty(compare_Edge_Edge_corr) && (range(orig_Edge_Edge_corr) > 0.01)
    histogram(orig_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
	set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_Edge_Edge_corr) < 0.01)
    histogram(compare_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
	set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Edge_Edge_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% GM_fd corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against overall movement.
nexttile
hold on
title('FD to GM Correlation', 'Interpreter', 'none')
histogram(orig_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_FD_GM_corr)
    histogram(compare_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(denoised_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% GM_dvars corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against global signal changes (including large
% movement)
nexttile
hold on
title('DVARS to GM Correlation', 'Interpreter', 'none')
histogram(orig_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_DVARS_GM_corr)
    histogram(compare_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(denoised_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% GM_Outbrain corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('Outbrain Correlation', 'Interpreter', 'none')
if ~isempty(compare_Outbrain_Outbrain_corr) && (range(orig_Outbrain_Outbrain_corr) > 0.01)
    histogram(orig_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_Outbrain_Outbrain_corr) < 0.01)
    histogram(compare_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Outbrain_Outbrain_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% GM_WMCSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('WMCSF Correlation', 'Interpreter', 'none')
if ~isempty(compare_WMCSF_WMCSF_corr) && (range(orig_WMCSF_WMCSF_corr) > 0.01)
    histogram(orig_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_WMCSF_WMCSF_corr) < 0.01)
    histogram(compare_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_WMCSF_WMCSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off


% GM_CSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('CSF Correlation', 'Interpreter', 'none')
if ~isempty(compare_CSF_CSF_corr) && (range(orig_CSF_CSF_corr) > 0.01)
    histogram(orig_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_CSF_CSF_corr) < 0.01)
    histogram(compare_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_CSF_CSF_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% Suscept corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that. Want to see the denoised and orig be significantly different
% does not need to be as drastic as others, as susceptibility regions
% include GM and can be a bit subjective
nexttile
hold on
title('Susceptibility Correlation', 'Interpreter', 'none')
if ~isempty(compare_Suscept_Suscept_corr) && (range(orig_Suscept_Suscept_corr) > 0.01)
    histogram(orig_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_Suscept_Suscept_corr) < 0.01)
    histogram(compare_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Suscept_Suscept_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% GM_NotGM corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that. Want to see the denoised and orig be significantly different
nexttile
hold on
title('NotGM Correlation', 'Interpreter', 'none')
if ~isempty(compare_NotGM_NotGM_corr) && (range(orig_NotGM_NotGM_corr) > 0.01)
    histogram(orig_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_NotGM_NotGM_corr) < 0.01)
    histogram(compare_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_NotGM_NotGM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off


% GM_GM corr to orig: Want to see signal correlations maintained, in comparison to
% GM_NotGM corr, but still reduced (removed erroneous/global correlations).
% get color order right to match other plots
nexttile
hold on
title('GM Correlation', 'Interpreter', 'none')
if ~isempty(compare_GM_GM_corr) && (range(orig_GM_GM_corr) > 0.01)
    histogram(orig_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(compare_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
elseif (range(orig_GM_GM_corr) < 0.01)
    histogram(compare_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
else
    histogram(orig_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_GM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    lgd = legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
end
%xline(0, 'DisplayName','r = 0')
xlabel('Correlation (r)')
ylabel('Frequency')
xlim([-1, 1]);
hold off

% for single subject, you have one more plot to make
if ~isempty(denoised_GM_mean)
    % Plot GS changes before and after Denoising - it is easier to
    % detrend it all first for comparison
    nexttile([1 3])
    hold on
    title('Mean GM Signal (Detrended)', 'Interpreter', 'none')
    plot(detrend(orig_GM_mean,2), 'LineWidth', 1.5)
    plot(detrend(compare_GM_mean,2), 'LineWidth', 1.5)
    plot(detrend(denoised_GM_mean,2), 'LineWidth', 1.5)
    lgd = legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none');
    set(lgd, 'FontSize', lgd_font_size, 'Box', 'off');
    xlabel('Timepoints')
    ylabel('Signal')
    xlim([0, length(orig_GM_mean)])
    hold off
end

fprintf('Saving Figure\n')
save_figure_robust(fig, qc_plot_dest, 300)

%exportgraphics(t, qc_plot_dest, 'Resolution', 300)
end