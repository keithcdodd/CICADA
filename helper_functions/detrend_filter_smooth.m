function [cleaned_file] = detrend_filter_smooth(file, funcmask, output_dir, smoothing_kernel, fpass, detrend_degree)
% function to apply detrending (to 2nd polynomial by default), bandpass, and apply gaussian
% smoothing kernel, then write file to output dir
% fpass should be an array of two numbers, representing bounds of Hz to
% bandpass, e.g. [0.008,0.15]
% One could do low pass filtering, for example, with [0,0.15].

% Set up: read file and convert the requested physical smoothing kernel to
% a per-axis voxel-space sigma.
[~, file_name, ~] = fileparts(file);
file_orig_data = niftiread(file);
file_orig_data_info = niftiinfo(file);
tr = file_orig_data_info.PixelDimensions(4); % grab tr
N = file_orig_data_info.ImageSize(4); % number of samples
fs = 1/tr; % grab sampling rate for potential bandpass
funcmask_data = niftiread(funcmask);
    
voxel_sizes_mm = file_orig_data_info.PixelDimensions(1:3);
[sigma_vox, sigma_mm, voxel_sizes_mm] = ...
    cicada_fwhm_mm_to_sigma_vox(smoothing_kernel, voxel_sizes_mm);
fprintf(['  Output smoothing geometry: requested FWHM = %.12g mm; ', ...
    'voxel sizes = [%.12g %.12g %.12g] mm; sigma = ', ...
    '%.12g mm = [%.12g %.12g %.12g] voxels.\n'], ...
    smoothing_kernel, voxel_sizes_mm, sigma_mm, sigma_vox);

if ~isequal(size(funcmask_data), size(file_orig_data(:,:,:,1)))
    fprintf('   Funcmask size does not match Data size...\n')
    return
end


%% Temporal processing

% Preserve historical behavior when optional arguments are absent
if exist('fpass', 'var') ~= 1
    fpass = [];
end

if exist('detrend_degree', 'var') ~= 1
    detrend_degree = [];
end

file_orig_data_2D = reshape( ...
    file_orig_data, [], size(file_orig_data,4));

% Preserve historical CICADA arithmetic/order exactly.
% Calculate voxel means while data are [nVoxels x nTime].
file_means_2D = mean(file_orig_data_2D, 2);

[data_temporal, ~, temporal_info] = ...
    cicada_temporal_transform( ...
        file_orig_data_2D', ...
        tr, ...
        fpass, ...
        detrend_degree, ...
        file_means_2D');

% Historical BOLD mean restoration remains outside the
% mathematical temporal transform T.
filtered_signal_2D = data_temporal' + file_means_2D;

detrended = temporal_info.detrended;
bp = temporal_info.filtered;

filtered_signal = reshape( ...
    filtered_signal_2D, ...
    size(file_orig_data, 1), ...
    size(file_orig_data, 2), ...
    size(file_orig_data, 3), ...
    size(file_orig_data, 4));

% smooth file and rename appropriately
cleaned_file = cicada_spatial_smooth_write( ...
    filtered_signal, ...
    funcmask_data, ...
    file_orig_data_info, ...
    file_name, ...
    output_dir, ...
    smoothing_kernel, ...
    sigma_vox, ...
    bp, ...
    detrended);

end
