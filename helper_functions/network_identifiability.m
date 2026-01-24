function ident_table = network_identifiability(template_file, cleaned_file, compare_file, orig_file, funcmask, GM_mask, smooth_kern, network_names, output_dir)
%NETWORK_IDENTIFIABILITY
% Quantifies functional network identifiability and generates visualization
% maps for comparison of denoising pipelines.
%
% For each template-defined network (seed parcel), a representative time
% series is extracted (mean or PC1) and correlated with all voxels within
% the functional brain mask. Correlation maps are Fisher z-transformed and
% assembled into 4D network-specific z-maps.
%
% Quantitative identifiability metrics are computed using a joint criterion:
%   (1) Evidence threshold: strongest network z-value > z_thr
%   (2) Separability threshold: (strongest – second-strongest) z-value > delta
%
% Metrics are reported as voxel fractions within:
%   • Global gray matter (GM ∩ functional mask)
%   • Each network’s seed parcel
%
% For each region, the following are computed:
%   • PassThresh: fraction of voxels exceeding the evidence threshold
%   • PassDeltaGivenPassThresh: fraction exceeding the separability threshold
%     conditional on passing the evidence threshold
%   • PassBoth: fraction satisfying both criteria simultaneously
%
% Results are written to a human-readable CSV with rows corresponding to
% global GM and individual networks, and columns corresponding to cleaned,
% comparison, and original data, as well as their pairwise differences.
%
% In addition, thresholded dominant-network label maps are generated for
% visualization and written as a multi-volume NIfTI file (cleaned, compare,
% original).
%
% Notes:
% - Quantitative metrics are computed within GM ∩ functional mask.
% - Visualization maps are thresholded independently for display purposes.
% - Networks are indexed internally as iNet = 1..K to support arbitrary
%   template label values.


    % ----------------------------
    % User-tunable discriminability params, but these settings generally
    % work well
    % ----------------------------
    use_pc1 = false;      % false => mean time series (default). true => PC1 from seed voxels.
    z_thr   = 0.6;       % optional: require winner/correct evidence above this (Fisher-z); set 0 to disable
    delta   = 0.15;       % margin threshold for "high confidence" (Fisher-z units)
    % The idea above is we want to balance/reward both when more voxels are
    % above a threshold for matching a network AND that they can be
    % differentiated from other networks reasonably (moderate instead of large delta given
    % networks are not entirely independent in practice)

    % Ensure output directory exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % get template dir
    [template_dir,~,~] = fileparts(template_file);

    % Temporary resampled template path
    resampled_template = fullfile(template_dir, 'resampled_template.nii.gz');

    % resample template to functional space
    call_fsl(['flirt -in ', template_file,' -ref ', cleaned_file,' -applyxfm -usesqform -interp nearestneighbour -out ', resampled_template]);

    % Load volumes
    template_vol = niftiread(resampled_template);
    funcmask_vol = niftiread(funcmask);
    GM_vol       = niftiread(GM_mask);
    cleaned_vol  = niftiread(cleaned_file);
    compare_vol  = niftiread(compare_file);
    orig_vol     = niftiread(orig_file);

    % Dimension checks
    assert(all(size(cleaned_vol) == size(orig_vol)),   'Cleaned and original dimensions must match');
    assert(all(size(cleaned_vol) == size(compare_vol)), 'Cleaned and compare dimensions must match');
    assert(all(size(cleaned_vol(:,:,:,1)) == size(funcmask_vol(:,:,:,1))), 'Cleaned and mask dimensions must match');
    assert(all(size(template_vol(:,:,:,1)) == size(cleaned_vol(:,:,:,1))), 'Template and fMRI must be in same space');

    % ----------------------------
    % Smooth cleaned/compare/orig
    % ----------------------------
    cleaned_vol_smoothed = zeros(size(cleaned_vol), 'like', cleaned_vol);
    compare_vol_smoothed = zeros(size(compare_vol), 'like', compare_vol);
    orig_vol_smoothed    = zeros(size(orig_vol),   'like', orig_vol);

    cleaned_vol_info = niftiinfo(cleaned_file);
    voxel_size = mean(cleaned_vol_info.PixelDimensions(1:3));

    fwhm_mm     = smooth_kern;
    fwhm_voxels = fwhm_mm / voxel_size;
    sigma       = fwhm_voxels / 2.355;

    GM_eval = logical(single(GM_vol > 0) .* single(funcmask_vol > 0));


    for t = 1:size(cleaned_vol, 4)
        cleaned_vol_smoothed(:,:,:,t) = imgaussfilt3(cleaned_vol(:,:,:,t), sigma);
        compare_vol_smoothed(:,:,:,t) = imgaussfilt3(compare_vol(:,:,:,t), sigma);
        orig_vol_smoothed(:,:,:,t)    = imgaussfilt3(orig_vol(:,:,:,t), sigma);
    end

    cleaned_vol = cleaned_vol_smoothed;
    compare_vol = compare_vol_smoothed;
    orig_vol    = orig_vol_smoothed;

    % ----------------------------
    % Evaluation mask: GM ∩ funcmask
    % ----------------------------
    % Masks
    vis_mask  = (funcmask_vol > 0);   % for building/visualizing z maps
    eval_mask = GM_eval;


    % Mask template by eval mask so seeds are only considered where you evaluate
    % template_vol(~eval_mask) = 0;

    % Network IDs (seed labels)
    network_ids = unique(template_vol(template_vol > 0 & vis_mask));
    network_ids = network_ids(:)';
    K = numel(network_ids);

    % Basic sizes
    [nx, ny, nz, nt] = size(cleaned_vol);

    % Preallocate 4D z-maps with robust indexing (iNet = 1..K)
    z_map_4D_cleaned = nan(nx, ny, nz, K, 'single');
    z_map_4D_compare = nan(nx, ny, nz, K, 'single');
    z_map_4D_orig    = nan(nx, ny, nz, K, 'single');

    % Helper to reshape 4D data into [voxels x time] given a mask
    vol_2d = @(vol, msk) reshape(vol(repmat(msk, 1, 1, 1, nt)), [], nt);


    % ----------------------------
    % Iterate over networks (seed parcels)
    % ----------------------------
    for iNet = 1:K
        net_id = network_ids(iNet);

        % Seed mask for this network (within eval mask)
        seed_mask_ts = (template_vol == net_id) & vis_mask;   % for seed time series

        % Extract voxel time series within current seed
        cleaned_net_ts = vol_2d(cleaned_vol, seed_mask_ts);
        compare_net_ts = vol_2d(compare_vol, seed_mask_ts);
        orig_net_ts    = vol_2d(orig_vol, seed_mask_ts);

        % Representative time series (mean or PC1)
        if use_pc1
            % PCA expects observations in rows, variables in columns; we want time as observations
            [~, score_cleaned] = pca(cleaned_net_ts', 'NumComponents', 1);
            [~, score_compare] = pca(compare_net_ts', 'NumComponents', 1);
            [~, score_orig]    = pca(orig_net_ts',    'NumComponents', 1);

            pc_cleaned = score_cleaned(:,1);
            pc_compare = score_compare(:,1);
            pc_orig    = score_orig(:,1);
        else
            pc_cleaned = mean(cleaned_net_ts, 1)';
            pc_compare = mean(compare_net_ts, 1)';
            pc_orig    = mean(orig_net_ts, 1)';
        end

        % Correlate with all visual mask (funcmask)
        cleaned_all_ts = vol_2d(cleaned_vol, vis_mask);
        compare_all_ts = vol_2d(compare_vol, vis_mask);
        orig_all_ts    = vol_2d(orig_vol,    vis_mask);

        corr_cleaned = corr(pc_cleaned, cleaned_all_ts', 'Rows', 'pairwise');
        corr_compare = corr(pc_compare, compare_all_ts', 'Rows', 'pairwise');
        corr_orig    = corr(pc_orig,    orig_all_ts', 'Rows', 'pairwise');

        % Fisher Z
        z_cleaned = atanh(corr_cleaned);
        z_compare = atanh(corr_compare);
        z_orig    = atanh(corr_orig);

        % Reconstruct full-brain z maps
        z_map_cleaned = nan(nx, ny, nz);
        z_map_compare = nan(nx, ny, nz);
        z_map_orig    = nan(nx, ny, nz);

        idx = find(vis_mask);
        z_map_cleaned(idx) = z_cleaned;
        z_map_compare(idx) = z_compare;
        z_map_orig(idx)    = z_orig;

        % Handle inf
        z_map_cleaned(isinf(z_map_cleaned)) = NaN;
        z_map_compare(isinf(z_map_compare)) = NaN;
        z_map_orig(isinf(z_map_orig))       = NaN;

        % Store in 4D stacks (indexed by iNet, not net_id)
        z_map_4D_cleaned(:,:,:,iNet) = single(z_map_cleaned);
        z_map_4D_compare(:,:,:,iNet) = single(z_map_compare);
        z_map_4D_orig(:,:,:,iNet)    = single(z_map_orig);

        % ----------------------------
        % Includ per-network discriminability (within this seed parcel):
        % "Correct-vs-best-other" margin: z_correct - max(z_other_networks)
        % computed later after all networks are built (needs full z_map_4D_*).
        % For now, store placeholders; we'll fill after the loop.
        % ----------------------------
        mean_margin_seed_cleaned = NaN; frac_highconf_seed_cleaned = NaN;
        mean_margin_seed_compare = NaN; frac_highconf_seed_compare = NaN;
        mean_margin_seed_orig    = NaN; frac_highconf_seed_orig    = NaN;

        % Network name lookup (robust)
        if nargin >= 8 && ~isempty(network_names) && net_id <= numel(network_names) && ~isempty(network_names{net_id})
            net_name = network_names{net_id};
        else
            net_name = sprintf('Net_%d', net_id);
        end
    end


    % testing diagnostics, can delete later
    % ----------------------------
    % Diagnostic: threshold vs margin behavior
    % ----------------------------
    fprintf('\n--- Network Identifiability Diagnostics ---\n');
    fprintf('z_thr = %.2f, delta = %.2f\n', z_thr, delta);
    
    Zc_tmp = z_map_4D_cleaned; Zc_tmp(~isfinite(Zc_tmp)) = -Inf;
    Zp_tmp = z_map_4D_compare; Zp_tmp(~isfinite(Zp_tmp)) = -Inf;
    Zo_tmp = z_map_4D_orig;    Zo_tmp(~isfinite(Zo_tmp)) = -Inf;
    
    labels = {'Cleaned', 'Compare', 'Orig'};
    Zall = {Zc_tmp, Zp_tmp, Zo_tmp};
    
    for ii = 1:3
        Z = Zall{ii};
    
        % Top1 and top2
        Zs = sort(Z, 4, 'descend');
        top1 = Zs(:,:,:,1);
        top2 = Zs(:,:,:,2);
    
        % Apply eval mask
        t1 = top1(eval_mask);
        t2 = top2(eval_mask);
    
        % Coverage at z_thr
        frac_above = mean(t1 > z_thr, 'omitnan');
    
        % Conditional separability
        if any(t1 > z_thr)
            frac_sep = mean((t1(t1 > z_thr) - t2(t1 > z_thr)) > delta, 'omitnan');
        else
            frac_sep = NaN;
        end

        % combined
        frac_both = mean( (t1 > z_thr) & ((t1 - t2) > delta), 'omitnan' );

    
        fprintf('%s:\n', labels{ii});
        fprintf('  P(top1 > z_thr)           = %.4f\n', frac_above);
        fprintf('  P(margin > delta | pass) = %.4f\n', frac_sep);
        fprintf('  P(pass both)             = %.4f\n', frac_both);
    end
    fprintf('------------------------------------------\n\n');

    % ----------------------------
    % Quantitative metrics for CSV (readable format)
    % ----------------------------
    
    % Ensure masks are logical
    vis_mask  = logical(vis_mask);
    eval_mask = logical(eval_mask);
    
    % Compute global metrics (GM expanded mask) for each pipeline
    [gm_thr_cleaned, gm_dlt_cleaned, gm_both_cleaned] = compute_global_pass_metrics(z_map_4D_cleaned, eval_mask, z_thr, delta);
    [gm_thr_compare, gm_dlt_compare, gm_both_compare] = compute_global_pass_metrics(z_map_4D_compare, eval_mask, z_thr, delta);
    [gm_thr_orig,    gm_dlt_orig,    gm_both_orig]    = compute_global_pass_metrics(z_map_4D_orig,    eval_mask, z_thr, delta);
    
    % Compute per-network seed metrics for each pipeline
    [seed_thr_cleaned, seed_dlt_cleaned, seed_both_cleaned] = compute_seed_pass_metrics(z_map_4D_cleaned, template_vol, eval_mask, network_ids, z_thr, delta);
    [seed_thr_compare, seed_dlt_compare, seed_both_compare] = compute_seed_pass_metrics(z_map_4D_compare, template_vol, eval_mask, network_ids, z_thr, delta);
    [seed_thr_orig,    seed_dlt_orig,    seed_both_orig]    = compute_seed_pass_metrics(z_map_4D_orig,    template_vol, eval_mask, network_ids, z_thr, delta);
    
    % Build row names (GM first, then each network)
    row_names = {};
    row_values = [];  % Nx5: [Cleaned, Compare, Orig, C-Comp, C-Orig]
    
    % GM rows
    row_names(end+1,1) = {"GM Overall Score"};
    row_values(end+1,:) = [gm_both_cleaned, gm_both_compare, gm_both_orig, gm_both_cleaned-gm_both_compare, gm_both_cleaned-gm_both_orig];
    
    row_names(end+1,1) = {"GM PassThresh"};
    row_values(end+1,:) = [gm_thr_cleaned, gm_thr_compare, gm_thr_orig, gm_thr_cleaned-gm_thr_compare, gm_thr_cleaned-gm_thr_orig];
    
    row_names(end+1,1) = {"GM PassDeltaGivenPassThresh"};
    row_values(end+1,:) = [gm_dlt_cleaned, gm_dlt_compare, gm_dlt_orig, gm_dlt_cleaned-gm_dlt_compare, gm_dlt_cleaned-gm_dlt_orig];
    
    % Per-network rows
    for iNet = 1:numel(network_ids)
        net_id = network_ids(iNet);
    
        % Network display name
        if nargin >= 8 && ~isempty(network_names) && net_id <= numel(network_names) && ~isempty(network_names{net_id})
            net_name = network_names{net_id};
        else
            net_name = sprintf('Net_%d', net_id);
        end
    
        % PassBoth, which is the overall score
        row_names(end+1,1) = {sprintf('%s Overall Score', net_name)};
        row_values(end+1,:) = [seed_both_cleaned(iNet), seed_both_compare(iNet), seed_both_orig(iNet), ...
                               seed_both_cleaned(iNet)-seed_both_compare(iNet), seed_both_cleaned(iNet)-seed_both_orig(iNet)];
    
        % PassThresh
        row_names(end+1,1) = {sprintf('%s PassThresh', net_name)};
        row_values(end+1,:) = [seed_thr_cleaned(iNet), seed_thr_compare(iNet), seed_thr_orig(iNet), ...
                               seed_thr_cleaned(iNet)-seed_thr_compare(iNet), seed_thr_cleaned(iNet)-seed_thr_orig(iNet)];
    
        % PassDeltaGivenPassThresh
        row_names(end+1,1) = {sprintf('%s PassDeltaGivenPassThresh', net_name)};
        row_values(end+1,:) = [seed_dlt_cleaned(iNet), seed_dlt_compare(iNet), seed_dlt_orig(iNet), ...
                               seed_dlt_cleaned(iNet)-seed_dlt_compare(iNet), seed_dlt_cleaned(iNet)-seed_dlt_orig(iNet)];
    end
    
    % Create final readable table
    T = table(row_names, ...
              row_values(:,1), row_values(:,2), row_values(:,3), row_values(:,4), row_values(:,5), ...
              'VariableNames', {'Metric', 'Cleaned', 'Compare', 'Orig', 'CleanedMinusCompare', 'CleanedMinusOrig'});
    
    % Write CSV (replaces old format)
    writetable(T, fullfile(output_dir, 'network_identifiability_results.csv'));
    
    % Return table
    ident_table = T;

    % ----------------------------
    % Visualization outputs
    % Dominant network label maps from thresholded z maps
    % ----------------------------

    % Threshold for visualization  
    viz_z = 0.67; % picking something strict for easier visualization

    % Build thresholded 4D arrays for max() labeling (replace NaN/Inf with 0)
    Zc = z_map_4D_cleaned; Zc(~isfinite(Zc)) = 0; Zc(Zc < viz_z) = 0;
    Zp = z_map_4D_compare; Zp(~isfinite(Zp)) = 0; Zp(Zp < viz_z) = 0;
    Zo = z_map_4D_orig;    Zo(~isfinite(Zo)) = 0; Zo(Zo < viz_z) = 0;

    % Argmax across networks gives index 1..K; convert to label values network_ids(index)
    [~, maxIdx_cleaned] = max(Zc, [], 4);
    [~, maxIdx_compare] = max(Zp, [], 4);
    [~, maxIdx_orig]    = max(Zo, [], 4);

    zeroMask_cleaned = all(Zc == 0, 4);
    zeroMask_compare = all(Zp == 0, 4);
    zeroMask_orig    = all(Zo == 0, 4);

    maxIdx_cleaned(zeroMask_cleaned) = 0;
    maxIdx_compare(zeroMask_compare) = 0;
    maxIdx_orig(zeroMask_orig)       = 0;

    % Convert iNet indices to actual network label IDs
    dominantNetworkMap_cleaned = zeros(nx, ny, nz, 'single');
    dominantNetworkMap_compare = zeros(nx, ny, nz, 'single');
    dominantNetworkMap_orig    = zeros(nx, ny, nz, 'single');

    dominantNetworkMap_cleaned(maxIdx_cleaned > 0) = network_ids(maxIdx_cleaned(maxIdx_cleaned > 0));
    dominantNetworkMap_compare(maxIdx_compare > 0) = network_ids(maxIdx_compare(maxIdx_compare > 0));
    dominantNetworkMap_orig(maxIdx_orig > 0)       = network_ids(maxIdx_orig(maxIdx_orig > 0));

    dominantNetworkMap_cleaned(~vis_mask) = 0;
    dominantNetworkMap_compare(~vis_mask) = 0;
    dominantNetworkMap_orig(~vis_mask)    = 0;

    dominantNetworkMaps = cat(4, dominantNetworkMap_cleaned, dominantNetworkMap_compare, dominantNetworkMap_orig);

    % Write nifti (3 volumes: cleaned, compare, orig)
    nifti_info = niftiinfo(funcmask);
    nifti_info.Datatype = 'single';
    nifti_info.PixelDimensions = [nifti_info.PixelDimensions, 1];
    nifti_info.ImageSize = [nifti_info.ImageSize, 3];
    nifti_basename = [output_dir, '/network_identifiability'];
    nifti_info.Filename = [output_dir, '/network_identifiability.nii.gz'];
    niftiwrite(single(dominantNetworkMaps), nifti_basename, nifti_info, 'Compressed', true);

end

% =========================================================================
% Helper functions
% =========================================================================

function [pass_thresh, pass_delta_given, pass_both] = compute_global_pass_metrics(Z4D, eval_mask, z_thr, delta)
    K = size(Z4D, 4);
    if K < 2
        pass_thresh = NaN; pass_delta_given = NaN; pass_both = NaN;
        return;
    end

    % Extract eval voxels into [Nvox x K]
    Zv = reshape(Z4D(repmat(eval_mask, 1, 1, 1, K)), [], K);
    Zv(~isfinite(Zv)) = NaN;

    % Rank networks per voxel (NaN -> -Inf for ranking)
    Zr = Zv;
    Zr(isnan(Zr)) = -Inf;
    Zs = sort(Zr, 2, 'descend');

    top1 = Zs(:,1);
    top2 = Zs(:,2);
    margin = top1 - top2;

    pass1 = (top1 > z_thr) & isfinite(margin);
    pass2 = pass1 & (margin > delta);

    pass_thresh = mean(pass1, 'omitnan');
    pass_both   = mean(pass2, 'omitnan');

    if any(pass1)
        pass_delta_given = mean(margin(pass1) > delta, 'omitnan');
    else
        pass_delta_given = NaN;
    end
end

function [pass_thresh, pass_delta_given, pass_both] = compute_seed_pass_metrics(Z4D, template_vol, eval_mask, network_ids, z_thr, delta)
    K = size(Z4D, 4);
    N = numel(network_ids);

    pass_thresh = nan(N,1);
    pass_delta_given = nan(N,1);
    pass_both = nan(N,1);

    if K < 2
        return;
    end

    Ztmp = Z4D;
    Ztmp(~isfinite(Ztmp)) = -Inf;

    for iNet = 1:N
        net_id = network_ids(iNet);
        m = (template_vol == net_id) & eval_mask;
        if ~any(m(:))
            continue;
        end

        zc_vol = Ztmp(:,:,:,iNet);
        zc = zc_vol(m);

        otherIdx = setdiff(1:K, iNet);
        z_other_max = -Inf(size(zc));
        for j = 1:numel(otherIdx)
            zo_vol = Ztmp(:,:,:,otherIdx(j));
            z_other_max = max(z_other_max, zo_vol(m));
        end

        margin = zc - z_other_max;

        pass1 = (zc > z_thr) & isfinite(margin);
        pass2 = pass1 & (margin > delta);

        pass_thresh(iNet) = mean(pass1, 'omitnan');
        pass_both(iNet)   = mean(pass2, 'omitnan');

        if any(pass1)
            pass_delta_given(iNet) = mean(margin(pass1) > delta, 'omitnan');
        else
            pass_delta_given(iNet) = NaN;
        end
    end
end
