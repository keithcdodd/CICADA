function ident_table = network_identifiability(template_file, cleaned_file, compare_file, orig_file, funcmask, GM_mask, smooth_kern, network_names, output_dir)
    % Ensure output directory exists
    if ~exist(output_dir, 'dir')
        mkdir(output_dir);
    end

    % get template dir:
    [template_dir,~,~] = fileparts(template_file);

    % Temporary resampled template path
    resampled_template = fullfile(template_dir, 'resampled_template.nii.gz');
    
    % fsl flirt to resample template if needed
    call_fsl(['flirt -in ', template_file,' -ref ', cleaned_file,' -applyxfm -usesqform -interp nearestneighbour -out ', resampled_template]);
    
    
    % Load volumes
    template_vol = niftiread(resampled_template);
    funcmask_vol = niftiread(funcmask);
    GM_vol = niftiread(GM_mask);
    cleaned_vol = niftiread(cleaned_file);
    compare_vol = niftiread(compare_file);
    orig_vol = niftiread(orig_file);


    % Get volume size and check dimensions
    assert(all(size(cleaned_vol) == size(orig_vol)), 'Cleaned and original dimensions must match');
    assert(all(size(cleaned_vol) == size(compare_vol)), 'Cleaned and compare dimensions must match');
    assert(all(size(cleaned_vol(:,:,:,1)) == size(funcmask_vol(:,:,:,1))), 'Cleaned and mask dimensions must match');
    assert(all(size(template_vol(:,:,:,1)) == size(cleaned_vol(:,:,:,1))), 'Template and fMRI must be in same space');

    % Smoothing cleaned and orig
    %fprintf('Smoothing functional data with kernel = %.2f mm\n', smooth_kern);
    cleaned_vol_smoothed = zeros(size(cleaned_vol), 'like', cleaned_vol);
    compare_vol_smoothed = zeros(size(compare_vol), 'like', compare_vol);
    orig_vol_smoothed = zeros(size(orig_vol), 'like', orig_vol);

    cleaned_vol_info = niftiinfo(cleaned_file);
    voxel_size = mean(cleaned_vol_info.PixelDimensions(1:3)); % to do the mm you want
    
    % sigma is about FWHMx / 2.355
    fwhm_mm = smooth_kern; % FWHMx
    fwhm_voxels = fwhm_mm / voxel_size;   
    sigma = fwhm_voxels / 2.355;          

    for t = 1:size(cleaned_vol, 4)
        cleaned_vol_smoothed(:,:,:,t) = imgaussfilt3(cleaned_vol(:,:,:,t), sigma);
        compare_vol_smoothed(:,:,:,t) = imgaussfilt3(compare_vol(:,:,:,t), sigma);
        orig_vol_smoothed(:,:,:,t) = imgaussfilt3(orig_vol(:,:,:,t), sigma);
    end

    cleaned_vol = cleaned_vol_smoothed;
    compare_vol = compare_vol_smoothed;
    orig_vol = orig_vol_smoothed;

    % mask template_vol by functional mask, for best results, look at data
    % within the GM mask
    orig_template_vol = template_vol;
    template_vol(funcmask_vol == 0) = 0;
    template_vol(GM_vol == 0) = 0;

    % Find network IDs
    network_ids = unique(template_vol(template_vol > 0));

    % Preallocate results
    results = [];

    % Iterate over networks
    for net_id = network_ids'
        mask = (orig_template_vol == net_id & funcmask_vol > 0); % within original network seed, in funcmask
        %other_mask = (template_vol > 0 & template_vol ~= net_id);
        other_mask = (funcmask_vol > 0 & orig_template_vol == 0); % outside of original network seeds, but within funcmask

        % Reshape 4D data into 2D: [voxels x time]
        [nx, ny, nz, nt] = size(cleaned_vol);
        vol_2d = @(vol, msk) reshape(vol(repmat(msk, 1, 1, 1, nt)), [], nt);

        % Get signals within current network
        cleaned_net_ts = vol_2d(cleaned_vol, mask);
        compare_net_ts = vol_2d(compare_vol, mask);
        orig_net_ts = vol_2d(orig_vol, mask);

        % Compute first PC
        % [~, score_cleaned, ~] = pca(cleaned_net_ts');
        % [~, score_compare, ~] = pca(compare_net_ts');
        % [~, score_orig, ~] = pca(orig_net_ts');
        % pc_cleaned = score_cleaned(:,1);
        % pc_compare = score_compare(:,1);
        % pc_orig = score_orig(:,1);

        % mean signal
        pc_cleaned = mean(cleaned_net_ts)';
        pc_compare = mean(compare_net_ts)';
        pc_orig = mean(orig_net_ts)';

        % Correlate PC with all network/funcmask voxels
        %all_net_mask = (template_vol > 0);
        all_net_mask = (funcmask_vol > 0);
        cleaned_all_ts = vol_2d(cleaned_vol, all_net_mask);
        compare_all_ts = vol_2d(compare_vol, all_net_mask);
        orig_all_ts = vol_2d(orig_vol, all_net_mask);

        % Correlation
        corr_cleaned = corr(pc_cleaned, cleaned_all_ts');
        corr_compare = corr(pc_compare, compare_all_ts');
        corr_orig = corr(pc_orig, orig_all_ts');

        % Fisher Z-transform
        z_cleaned = atanh(corr_cleaned);
        z_compare = atanh(corr_compare);
        z_orig = atanh(corr_orig);

        % Reconstruct full-brain z-score maps
        z_map_cleaned = nan(nx, ny, nz);
        z_map_compare = nan(nx, ny, nz);
        z_map_orig = nan(nx, ny, nz);

        idx = find(all_net_mask);
        z_map_cleaned(idx) = z_cleaned;
        z_map_compare(idx) = z_compare;
        z_map_orig(idx) = z_orig;

        % convert NaN to 0, inf to 1
        z_map_cleaned(isnan(z_map_cleaned)) = 0;
        z_map_cleaned(isinf(z_map_cleaned)) = 10; % just not inf
        z_map_compare(isnan(z_map_compare)) = 0;
        z_map_compare(isinf(z_map_compare)) = 10; % just not inf
        z_map_orig(isnan(z_map_orig)) = 0;
        z_map_orig(isinf(z_map_orig)) = 10; % just not inf

        % Save z-score maps
        z_map_4D_cleaned(:,:,:,net_id) = z_map_cleaned;
        z_map_4D_compare(:,:,:,net_id) = z_map_compare;
        z_map_4D_orig(:,:,:,net_id) = z_map_orig;

        % absolute value
        abs_z_cleaned = abs(z_map_cleaned);
        abs_z_compare = abs(z_map_compare);
        abs_z_orig = abs(z_map_orig);

        % Compute mean absolute z-values within and outside network
        in_abs_z_cleaned = abs_z_cleaned(mask);
        out_abs_z_cleaned = abs_z_cleaned(other_mask);

        in_abs_z_compare = abs_z_compare(mask);
        out_abs_z_compare = abs_z_compare(other_mask);

        in_abs_z_orig = abs_z_orig(mask);
        out_abs_z_orig = abs_z_orig(other_mask);
        
        % p-value 0.5 (50% chance) cut off
        cut_off = 0.67; %Z=0.67 is p=0.50 with two tails
        abs_z_cleaned_cut_off = abs_z_cleaned(abs_z_cleaned > cut_off);
        abs_z_compare_cut_off = abs_z_compare(abs_z_compare > cut_off);
        abs_z_orig_cut_off = abs_z_orig(abs_z_orig > cut_off);
        
        % Can use a cut off if desired
        cut_off = -0.01;
        mean_in_cleaned = mean(in_abs_z_cleaned(in_abs_z_cleaned > cut_off), 'omitnan');
        mean_out_cleaned = mean(out_abs_z_cleaned(out_abs_z_cleaned > cut_off), 'omitnan');
        ratio_cleaned = mean_in_cleaned / mean_out_cleaned;

        mean_in_compare = mean(in_abs_z_compare(in_abs_z_compare > cut_off), 'omitnan');
        mean_out_compare = mean(out_abs_z_compare(out_abs_z_compare > cut_off), 'omitnan');
        ratio_compare = mean_in_compare / mean_out_compare;

        mean_in_orig = mean(in_abs_z_orig(in_abs_z_orig > cut_off), 'omitnan');
        mean_out_orig = mean(out_abs_z_orig(out_abs_z_orig > cut_off), 'omitnan');
        ratio_orig = mean_in_orig / mean_out_orig;

        % Store results
        results = [results; table(net_id, network_names{net_id}, ratio_cleaned / ratio_compare, ratio_cleaned / ratio_orig, ...
                                  mean_in_cleaned, mean_out_cleaned, ratio_cleaned, ...
                                  mean_in_compare, mean_out_compare, ratio_compare, ...
                                  mean_in_orig, mean_out_orig, ratio_orig)];
    end

    % Rename columns
    results.Properties.VariableNames = {'NetworkID', 'Network_Names', 'Cleaned_Compare_Ratio', 'Cleaned_Orig_Ratio', ...
        'MeanIn_Cleaned', 'MeanOut_Cleaned', 'Ratio_Cleaned', ...
        'MeanIn_Compare', 'MeanOut_Compare', 'Ratio_Compare', ...
        'MeanIn_Orig', 'MeanOut_Orig', 'Ratio_Orig'};

    % Save table
    writetable(results, fullfile(output_dir, 'network_identifiability_results.csv'));

    % Want to make a 3D file that shows each network
    z_map_4D_cleaned_abs_cut_off = abs(z_map_4D_cleaned);
    z_map_4D_cleaned_abs_cut_off(z_map_4D_cleaned_abs_cut_off < 0.67) = 0; % z = 0.67 is 50% (p=0.50)

    z_map_4D_compare_abs_cut_off = abs(z_map_4D_compare);
    z_map_4D_compare_abs_cut_off(z_map_4D_compare_abs_cut_off < 0.67) = 0; % z = 0.67 is 50% (p=0.50)

    z_map_4D_orig_abs_cut_off = abs(z_map_4D_orig);
    z_map_4D_orig_abs_cut_off(z_map_4D_orig_abs_cut_off < 0.67) = 0; % z = 0.67 is 50% (p=0.50)

    % OK, now, find the max z-value across the 4th dimension, and label
    % voxel as that network index, if all 0, mark 0:
    % Get the maximum value and index across the 4th dimension
    [~, maxIdx_cleaned] = max(z_map_4D_cleaned_abs_cut_off, [], 4);
    [~, maxIdx_compare] = max(z_map_4D_compare_abs_cut_off, [], 4);
    [~, maxIdx_orig] = max(z_map_4D_orig_abs_cut_off, [], 4);
    
    % Create a logical mask where all values across the 4th dim are zero
    zeroMask_cleaned = all(z_map_4D_cleaned_abs_cut_off == 0, 4);
    zeroMask_compare = all(z_map_4D_compare_abs_cut_off == 0, 4);
    zeroMask_orig = all(z_map_4D_orig_abs_cut_off == 0, 4);
    
    % Assign zero to those voxels in the result where all networks were zero
    maxIdx_cleaned(zeroMask_cleaned) = 0;
    maxIdx_compare(zeroMask_compare) = 0;
    maxIdx_orig(zeroMask_orig) = 0;
    
    % The result is a 3D image where each voxel has the index of the most dominant network,
    % or 0 if no network was dominant
    dominantNetworkMap_cleaned = maxIdx_cleaned;
    dominantNetworkMap_compare = maxIdx_compare;
    dominantNetworkMap_orig = maxIdx_orig;

    dominantNetworkMaps = cat(4, dominantNetworkMap_cleaned, dominantNetworkMap_compare, dominantNetworkMap_orig);

    % write nifti
    nifti_info = niftiinfo(funcmask);
    nifti_info.Datatype = 'single';
    nifti_info.PixelDimensions = [nifti_info.PixelDimensions, 1];
    nifti_info.ImageSize = [nifti_info.ImageSize, 3]; % cleaned, compare, orig
    nifti_info.Filename = [output_dir, '/network_identifiability.nii.gz'];
    niftiwrite(single(dominantNetworkMaps), nifti_info.Filename, nifti_info, 'Compressed', true);

    % Return results
    ident_table = results;
end

