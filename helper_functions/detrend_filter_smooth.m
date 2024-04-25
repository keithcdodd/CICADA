function [cleaned_file] = detrend_filter_smooth(file, funcmask, output_dir, smoothing_kernel, fpass)
% function to apply detrending (to 2nd polynomial) and apply gaussian
% smoothing kernel, then write file to output dir
% fpass should be an array of two numbers, representing bounds of Hz to
% bandpass, e.g. [0.008,0.15]
% One could do low pass filtering, for example, with [0,0.15].

% Set up: read file, and get conversion between smoothing kernel from mm to
% voxel size
[~, file_name, ~] = fileparts(file);
file_orig_data = niftiread(file);
file_orig_data_info = niftiinfo(file);
tr = file_orig_data_info.PixelDimensions(4); % grab tr
N = file_orig_data_info.ImageSize(4); % number of samples
T = N * tr; % total time scanned
fs = 1/tr; % grab sampling rate for potential bandpass
funcmask_data = niftiread(funcmask);
mm_div = mean(file_orig_data_info.PixelDimensions(1:3)); % to convert smoothing kernel for imgaussfilt3D

if size(funcmask_data) ~= size(file_orig_data(:,:,:,1))
    fprintf('   Funcmask size does not match Data size...\n')
    return
end


% detrend
fprintf('   Detrending to 2nd polynomial...\n')
file_orig_data_2D = reshape(file_orig_data, [], size(file_orig_data,4)); % to apply detrend, filtering, etc.
file_means_2D = mean(file_orig_data_2D,2); % so that you can maintain the mean, if a program (like spm) wants to keep it.
file_detrended_2D = detrend(file_orig_data_2D',2)'; % mean is gone, back can be added back in later
filtered_signal_2D = file_detrended_2D;

bp = 0;
if (exist('fpass', 'var') == 1) && (isa(fpass, 'double') == 1) && ~isempty(fpass) && (length(fpass) == 2)
    bp = 1; % Mark bandpass filtering down
    % OK, do bandpass!
    fprintf(['  Bandpassing Filtering at ', num2str(fpass(1)), ' ', num2str(fpass(2)), ' Hz...\n'])
    % bandpass, can use fft and ifft
    N = T/tr; F = 1/tr;
    df = 1/T; N_freq = N/2 + 1;
    f = F*(0:floor(N/2))/N; % to chart what frequencies we are at
    lower_phys_cutoff = round(0.008 / df) + 1; % start at freq 0 at position 1
    higher_phys_cutoff = round(0.15 / df) + 1;
    
    % create for loop to loop through each thing
    filtered_signal_2D = zeros(size(file_detrended_2D));
    for idx1 = 1:size(file_detrended_2D,1)
        freq_signal = fftshift(fft(file_detrended_2D(idx1,:))); % center point at floor(N/2+1) is 0 Hz
        low_freq_cutoffs = (floor(N/2)+1 - lower_phys_cutoff): (floor(N/2)+1 + lower_phys_cutoff);
        high_freq_cutoffs_neg = 1:(floor(N/2)+1 - (higher_phys_cutoff+1));
        high_freq_cutoffs_pos = (floor(N/2)+1 + (higher_phys_cutoff+1)):(N);
        freq_signal_filtered = freq_signal;
        freq_signal_filtered(low_freq_cutoffs) = 0;
        freq_signal_filtered(high_freq_cutoffs_neg) = 0;
        freq_signal_filtered(high_freq_cutoffs_pos) = 0;
        freq_signal_filtered = fftshift(freq_signal_filtered);
        filtered_signal = ifft(freq_signal_filtered);
        filtered_signal_2D(idx1, :) = filtered_signal;
    end
else
    fprintf('   Not Applying Bandpass filtering!\n')
end

filtered_signal_2D = filtered_signal_2D + file_means_2D; % add mean back in, good for visualization purposes and makes spm happier

filtered_signal = reshape(filtered_signal_2D, size(file_orig_data, 1), size(file_orig_data, 2), size(file_orig_data, 3), size(file_orig_data, 4));


% mask filtered file so gauss padding is OK
file_filtered_data = zeros(size(file_orig_data));
for idx2 = 1:size(file_orig_data,4)
    curr_file_data = filtered_signal(:,:,:, idx2);
    curr_file_data(funcmask_data == 0) = 0; % mask it
    curr_file_data(curr_file_data < 0.01) = 0; % and get rid of negative numbers, if they exist
    file_filtered_data(:,:,:, idx2) = curr_file_data;
end


% smooth
if smoothing_kernel ~= 0
    fprintf(['  Smoothing at ', num2str(smoothing_kernel), 'mm gauss...\n'])
    file_data = zeros(size(file_orig_data));
    % gaussian smooth and remask
    for idx2 = 1:size(file_orig_data,4)
        curr_file_data = imgaussfilt3(file_filtered_data(:,:,:, idx2), round(smoothing_kernel ./ mm_div)); % convert kernel from mm to voxel
        curr_file_data(funcmask_data == 0) = 0; % mask it
        curr_file_data(curr_file_data < 0.01) = 0; % and get rid of negative numbers, if they exist
        file_data(:,:,:, idx2) = curr_file_data;
    end

    if bp == 1
        % write to data dir and relabel cleaned file:
        niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(smoothing_kernel), '_bp_', file_name], file_orig_data_info, "Compressed", true)
        cleaned_file = [output_dir, '/s', num2str(smoothing_kernel), '_bp_', file_name, '.gz']; % update cleaned_file
    else
        % write to data dir and relabel cleaned file:
        niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(smoothing_kernel), '_', file_name], file_orig_data_info, "Compressed", true)
        cleaned_file = [output_dir, '/s', num2str(smoothing_kernel), '_', file_name, '.gz']; % update cleaned_file
    end
    

else
    fprintf('   Not Smoothing! \n')
    if bp == 1
        file_data = file_filtered_data; % not smoothed, only detrended and masked if smoothing kernel was set at 0
        % write to data dir and relabel cleaned file:
        niftiwrite(cast(file_data, 'single'), [output_dir, '/bp_', file_name], file_orig_data_info, "Compressed", true)
        cleaned_file = [output_dir, '/bp_', file_name, '.gz']; % update cleaned_file
    else
        file_data = file_filtered_data; % not smoothed, only detrended and masked if smoothing kernel was set at 0
        % write to data dir and relabel cleaned file:
        niftiwrite(cast(file_data, 'single'), [output_dir, '/', file_name], file_orig_data_info, "Compressed", true)
        cleaned_file = [output_dir, '/', file_name, '.gz']; % update cleaned_file
    end
  
end

end