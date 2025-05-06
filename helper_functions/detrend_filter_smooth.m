function [cleaned_file] = detrend_filter_smooth(file, funcmask, output_dir, smoothing_kernel, fpass, detrend_degree)
% function to apply detrending (to 2nd polynomial by default), bandpass, and apply gaussian
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
    
% sigma is about FWHMx / 2.355
voxel_size = mean(file_orig_data_info.PixelDimensions(1:3)); % to convert smoothing kernel for imgaussfilt3D
fwhm_mm = smoothing_kernel; % FWHMx
fwhm_voxels = fwhm_mm / voxel_size;   
sigma = fwhm_voxels / 2.355; 

if size(funcmask_data) ~= size(file_orig_data(:,:,:,1))
    fprintf('   Funcmask size does not match Data size...\n')
    return
end


% detrend, if asked to
detrended = 0; % whether or not we detrended
if (exist('detrend_degree', 'var') == 1) && (isa(detrend_degree, 'double') == 1) && ~isempty(detrend_degree)
    % detrend degree is provided
    if detrend_degree > 0
        % OK, now it is actually legit
        detrend_degree = round(detrend_degree); % round to nearest integer
        fprintf('   Detrending to %d polynomial...\n', detrend_degree)
        file_orig_data_2D = reshape(file_orig_data, [], size(file_orig_data,4)); % to apply detrend, filtering, etc.
        file_means_2D = mean(file_orig_data_2D,2); % so that you can maintain the mean, if a program (like spm) wants to keep it.
        file_detrended_2D = detrend(file_orig_data_2D',detrend_degree)'; % mean is gone, back can be added back in later
        filtered_signal_2D = file_detrended_2D;

        detrended = 1;

    else
        % if 0 or negative value, don't detrend
        fprintf('   Not Detrending!\n')

        % still remove the mean to be consistent
        file_orig_data_2D = reshape(file_orig_data, [], size(file_orig_data,4)); % to apply detrend, filtering, etc.
        file_means_2D = mean(file_orig_data_2D,2); % so that you can maintain the mean, if a program (like spm) wants to keep it.
        file_detrended_2D = file_means_2D ; % mean is gone, back can be added back in later. Not actually detrended, just keeping naming consistent.
        filtered_signal_2D = file_detrended_2D;

    end
else
    % if detrend degree is not given, or not a number, go to default of 2
    detrend_degree = 2;

    fprintf('   Detrending to %d polynomial...\n', detrend_degree)
    file_orig_data_2D = reshape(file_orig_data, [], size(file_orig_data,4)); % to apply detrend, filtering, etc.
    file_means_2D = mean(file_orig_data_2D,2); % so that you can maintain the mean, if a program (like spm) wants to keep it.
    file_detrended_2D = detrend(file_orig_data_2D',detrend_degree)'; % mean is gone, back can be added back in later
    filtered_signal_2D = file_detrended_2D;

    detrended = 1;
    
end



bp = 0;
if (exist('fpass', 'var') == 1) && (isa(fpass, 'double') == 1) && ~isempty(fpass) && (length(fpass) == 2)
    bp = 1; % Mark bandpass filtering down
    % OK, do bandpass!

    % bandpass, can use fft and ifft
    N = T/tr; F = 1/tr;
    df = 1/T; N_freq = N/2 + 1;
    f = F*(0:floor(N/2))/N; % to chart what frequencies we are at

    % Now set it up to check values
    low_hz = fpass(1);
    high_hz = fpass(2);

    if low_hz < 0
        fprintf('   Low Hz Cut Off for Bandpass Filtering is Below 0... setting to 0 Hz. \n')
        low_hz = 0;
    end

    if high_hz > f(end)
        fprintf(['   High Hz Cut Off For Bandpass Filtering is Above Highest Frequency... setting to highest frequency of ', num2str(f(end)) ,' Hz. \n'])
        high_hz = f(end);
    end

    if low_hz >= high_hz
        fprintf('   Error: Low Hz Cut Off is Equal to or Greater Than High Hz Cut Off. Not Applying Bandpass Filtering. \n')

    else
        % Now we can actually apply bp filtering
        fprintf(['  Bandpassing Filtering at ', num2str(low_hz), ' ', num2str(high_hz), ' Hz...\n'])
        lower_phys_cutoff = round(low_hz / df) + 1; % start at freq 0 at position 1
        higher_phys_cutoff = round(high_hz / df) + 1;
        
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
        curr_file_data = imgaussfilt3(file_filtered_data(:,:,:, idx2), sigma); % convert kernel from mm to voxel
        curr_file_data(funcmask_data == 0) = 0; % mask it
        curr_file_data(curr_file_data < 0.01) = 0; % and get rid of negative numbers, if they exist
        file_data(:,:,:, idx2) = curr_file_data;
    end

    if bp == 1
        if detrended == 1
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(smoothing_kernel), '_bp_d', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(smoothing_kernel), '_bp_d_', file_name, '.gz']; % update cleaned_file
        else
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(smoothing_kernel), '_bp_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(smoothing_kernel), '_bp_', file_name, '.gz']; % update cleaned_file
        end
    else
        if detrended == 1
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(smoothing_kernel), '_d_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(smoothing_kernel), '_d_', file_name, '.gz']; % update cleaned_file
        else
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(smoothing_kernel), '_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(smoothing_kernel), '_', file_name, '.gz']; % update cleaned_file
        end
        
    end
    

else
    fprintf('   Not Smoothing! \n')
    if bp == 1
        if detrended == 1
            file_data = file_filtered_data; % not smoothed, only detrended, bandpassed, and masked if smoothing kernel was set at 0
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/bp_d_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/bp_d_', file_name, '.gz']; % update cleaned_file
        else
            file_data = file_filtered_data; % not smoothed, only bandpassed and masked if smoothing kernel was set at 0
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/bp_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/bp_', file_name, '.gz']; % update cleaned_file
        end
    else
        if detrended == 1
            file_data = file_filtered_data; % not smoothed, only detrended and masked if smoothing kernel was set at 0
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/d_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/d_', file_name, '.gz']; % update cleaned_file
        else
            file_data = file_filtered_data; % not smoothed or detrended, just masked if smoothing kernel was set at 0
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/', file_name, '.gz']; % update cleaned_file
        end
        
    end
  
end

end