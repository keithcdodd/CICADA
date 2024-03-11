function plot_qc(denoised_Edge_GM_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, denoised_Outbrain_GM_corr, ...
    denoised_WMCSF_GM_corr, denoised_CSF_GM_corr, denoised_NotGM_GM_corr, denoised_GM_GM_autocorr, ...
    compare_Edge_GM_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, compare_Outbrain_GM_corr, ...
    compare_WMCSF_GM_corr, compare_CSF_GM_corr, compare_NotGM_GM_corr, compare_GM_GM_autocorr, ...
    orig_Edge_GM_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, orig_Outbrain_GM_corr, ...
    orig_WMCSF_GM_corr, orig_CSF_GM_corr, orig_NotGM_GM_corr, orig_GM_GM_autocorr, ...
    denoised_GM_mean, compare_GM_mean, orig_GM_mean, title_string, qc_plot_dest, denoised_name)

% There needs to be a decision on whether we are including GM mean signal
% (good for within subject comparison, worthless in group qc comparison)
% You can see if GM mean is empty for denoised

fprintf('Creating Figures\n')


if isempty(denoised_GM_mean)
    figure('Position', [50, 250, 2000, 900])
    t = tiledlayout(2,4, 'Padding', 'compact', 'TileSpacing', 'compact');
else
    figure('Position', [50, 50, 2000, 900])
    t = tiledlayout(3,3, 'Padding', 'compact', 'TileSpacing', 'compact');
end

title(t, title_string, 'Interpreter', 'none')

% GM_Edge corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('Edge-GM Corr ', 'Interpreter', 'none')
histogram(orig_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_Edge_GM_corr)
    histogram(compare_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_Edge_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off

% GM_fd corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against overall movement.
nexttile
hold on
title('FD-GM Corr ', 'Interpreter', 'none')
histogram(orig_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_FD_GM_corr)
    histogram(compare_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_FD_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off


% GM_dvars corr: Want to see a shift toward center for denoised
% (removal of noise/sudden GS signal shift in GM). Suggests it
% is robust against global signal changes (including large
% movement)
nexttile
hold on
title('DVARS-GM Corr ', 'Interpreter', 'none')
histogram(orig_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_DVARS_GM_corr)
    histogram(compare_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_DVARS_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off

% GM_Outbrain corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('Outbrain-GM Corr ', 'Interpreter', 'none')
histogram(orig_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_Outbrain_GM_corr)
    histogram(compare_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_Outbrain_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off

% GM_WMCSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('WMCSF-GM Corr ', 'Interpreter', 'none')
histogram(orig_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_WMCSF_GM_corr)
    histogram(compare_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_WMCSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off


% GM_CSF corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that
nexttile
hold on
title('CSF-GM Corr ', 'Interpreter', 'none')
histogram(orig_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_CSF_GM_corr)
    histogram(compare_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_CSF_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off

% GM_NotGM corr: Want to see a tighter distribution, centered near
% 0, if original is significantly right shifted, want to see less
% of that. Want to see the denoised and orig be significantly different
nexttile
hold on
title('NotGM-GM Corr ', 'Interpreter', 'none')
histogram(orig_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
if ~isempty(compare_NotGM_GM_corr)
    histogram(compare_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    histogram(denoised_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', 'compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(denoised_NotGM_GM_corr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'LineWidth', 1.5)
    legend('orig', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off


% GM_GM autocorr to orig: Want to see signal correlations maintained, in comparison to
% GM_NotGM corr, but still reduced (removed erroneous/global correlations).
% get color order right to match other plots
colorord = get(gca, 'colororder');

nexttile
hold on
title('GM-GM Auto Corr to Orig', 'Interpreter', 'none')
if isempty(compare_GM_GM_autocorr)
    histogram(denoised_GM_GM_autocorr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    legend(sprintf('%s', denoised_name), 'Interpreter', 'none')
else
    histogram(compare_GM_GM_autocorr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(2,:), 'LineWidth', 1.5)
    histogram(denoised_GM_GM_autocorr, 'Normalization', 'pdf', 'DisplayStyle', 'stairs', 'EdgeColor', colorord(3,:), 'LineWidth', 1.5)
    legend('compare', sprintf('%s', denoised_name), 'Interpreter', 'none')
end
hold off

% for single subject, you have one more plot to make
if ~isempty(denoised_GM_mean)
    % Plot GS changes before and after Denoising - it is easier to
    % detrend it all first for comparison
    nexttile
    hold on
    title('Mean GM Signal (Detrended)', 'Interpreter', 'none')
    plot(detrend(orig_GM_mean,2), 'LineWidth', 1.5)
    plot(detrend(compare_GM_mean,2), 'LineWidth', 1.5)
    plot(detrend(denoised_GM_mean,2), 'LineWidth', 1.5)
    legend('orig', 'compare', 'CICADA', 'Interpreter', 'none')
    xlim([0, length(orig_GM_mean)])
    hold off
end


fprintf('Saving Figure\n')
exportgraphics(t, qc_plot_dest, 'Resolution', 300)
end