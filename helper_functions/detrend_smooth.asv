function [cleaned_file] = detrend_smooth(file, funcmask, output_dir, smoothing_kernel)
% function to apply detrending (to 2nd polynomial) and apply gaussian
% smoothing kernel, then write file to output dir

% Set up: read file, and get conversion between smoothing kernel from mm to
% voxel size
[~, file_name, ~] = fileparts(file);
file_orig_data = niftiread(file);
file_orig_data_info = niftiinfo(file);
funcmask_data = niftiread(funcmask);
mm_div = mean(file_orig_data_info.PixelDimensions(1:3)); % to convert smoothing kernel for imgaussfilt3D

if size(funcmask_data) ~= size(file_orig_data(:,:,:,1))
    fprintf('Funcmask size does not match Data size...\n')
    return
end

% detrend
file_orig_data_2D = reshape(file_orig_data, [], size(file_orig_data,4)); % to apply detrend
file_means_2D = mean(file_orig_data_2D); % so that you can maintain the mean, if a program (like spm) wants to keep it.
file_detrended_2D = detrend(file_orig_data_2D,2) + file_means_2D; % detrend up to 2nd polynomial. Keeps mean, just in case, but linear and quadratic drift is removed
file_detrended = reshape(file_detrended_2D, size(file_orig_data, 1), size(file_orig_data, 2), size(file_orig_data, 3), size(file_orig_data, 4));

% mask file_detrended so gauss padding is OK
file_detrended_data = zeros(size(file_orig_data));
for idx2 = 1:size(file_orig_data,4)
    curr_file_data = file_detrended(:,:,:, idx2);
    curr_file_data(funcmask_data == 0) = 0; % mask it
    curr_file_data(curr_file_data < 0.01) = 0; % and get rid of negative numbers, if they exist
    file_detrended_data(:,:,:, idx2) = curr_file_data;
end
file_detrended = file_detrended_data; % now it should be masked by funcmask!

% smooth
if smoothing_kernel ~= 0
    file_data = zeros(size(file_orig_data));
    for idx2 = 1:size(file_orig_data,4)
        curr_file_data = imgaussfilt3(file_detrended(:,:,:, idx2), round(smoothing_kernel ./ mm_div)); % convert kernel from mm to voxel
        curr_file_data(funcmask_data == 0) = 0; % mask it
        curr_file_data(curr_file_data < 0.01) = 0; % and get rid of negative numbers, if they exist
        file_data(:,:,:, idx2) = curr_file_data;
    end
    % write to data dir and relabel cleaned file:
    niftiwrite(cast(file_data, 'single'), [output_dir, '/s_', file_name], file_orig_data_info, "Compressed", true)
    cleaned_file = [output_dir, '/s_', file_name, '.gz']; % update cleaned_file

else
    file_data = file_detrended; % not smoothed, only detrended and masked if smoothing kernel was set at 0
    % write to data dir and relabel cleaned file:
    niftiwrite(cast(file_data, 'single'), [data_dir, '/', file_name], file_orig_data_info, "Compressed", true)
    cleaned_file = [output_dir, '/', file_name, '.gz']; % update cleaned_file
end

% I cannot get masking to work in matlab?, so try applying it with call fsl
% instead, or just figure out why it is not working as expected...

end