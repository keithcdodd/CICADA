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
fs = 1/tr; % grab sampling rate for potential bandpass
funcmask_data = niftiread(funcmask);
    
% sigma is about FWHMx / 2.355
voxel_size = round(mean(file_orig_data_info.PixelDimensions(1:3))); % to convert smoothing kernel for imgaussfilt3D
fwhm_mm = smoothing_kernel; % FWHM
fwhm_voxels = fwhm_mm / voxel_size;   
sigma = fwhm_voxels / 2.355; 

if ~isequal(size(funcmask_data), size(file_orig_data(:,:,:,1)))
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
        file_detrended_2D = file_orig_data_2D - file_means_2D; % remove mean, but can be added back in later. Not actually detrended, just keeping naming consistent.
        filtered_signal_2D = file_detrended_2D;
        detrended = 0;

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
    bp = 0; % Will mark to 1 once filtering is actually applied later below
    % OK, do bandpass!

    % bandpass, can use butter
    nyq = fs/2;  % move this ABOVE so it exists here

    % Now set it up to check values
    low_hz = fpass(1);
    high_hz = fpass(2);

    if low_hz < 0
        fprintf('   Low Hz Cut Off for Bandpass Filtering is Below 0... setting to 0 Hz. \n')
        low_hz = 0;
    end

    if high_hz > nyq
        fprintf(['   High Hz Cut Off For Bandpass Filtering is Above Nyquist... setting to Nyquist of ', num2str(nyq) ,' Hz. \n'])
        high_hz = nyq;
    end


    % Determine which components are enabled (per your conventions)
    doHighPass = (low_hz > 0);      % enabled only if strictly above 0
    doLowPass  = (high_hz < nyq);   % enabled only if strictly below Nyquist

    
    % Only treat low>=high as an error if both HP and LP are enabled (true band-pass)
    if doHighPass && doLowPass && (low_hz >= high_hz)
        fprintf('   Error: low_hz >= high_hz while requesting a band-pass. Skipping filtering.\n');

    else
        % Now we can actually apply temporal filtering (Butterworth + zero-phase)
        filtered_signal_2D = file_detrended_2D;  % default: no further filtering
        
        % If neither high or low pass is enabled, skip filtering entirely
        if ~doHighPass && ~doLowPass
            fprintf('   Filter bounds imply no filtering (low<=0 and high>=Nyquist). Skipping.\n');
            % leave filtered_signal_2D as whatever it currently is (i.e., file_detrended_2D)
        else
            fprintf(['  Filtering at ', num2str(low_hz), ' ', num2str(high_hz), ' Hz...\n'])    
            filt_order = 2;  % common choice; increase to 4 for sharper roll-off (with caution)
            % for Butter filter!
        
            if ~doHighPass && doLowPass
                % Low-pass only
                Wn = high_hz / nyq;
                Wn = min(max(Wn, 1e-6), 1-1e-6); % clamping
                [b,a] = butter(filt_order, Wn, 'low');
        
            elseif doHighPass && ~doLowPass
                % High-pass only
                Wn = low_hz / nyq;
                Wn = min(max(Wn, 1e-6), 1-1e-6); % clamping
                [b,a] = butter(filt_order, Wn, 'high');
        
            else
                % Band-pass
                Wn = [low_hz high_hz] / nyq;
                Wn = min(max(Wn, 1e-6), 1-1e-6); % clamping
                [b,a] = butter(filt_order, Wn, 'bandpass');
            end
        
            % Mark that filtering truly occurred
            bp = 1;
        
            minN = 3 * (max(length(a), length(b)) - 1) + 1; % safeguarding an error for short timeseries
            if N <= minN
                error('Time series too short for filtfilt (N=%d, need > %d).', N, minN);
            end

            % Apply zero-phase filtering along the TIME dimension
            % file_detrended_2D is [nVox x nTime] -> transpose to [nTime x nVox] so filtfilt filters along time
            filtered_signal_2D = filtfilt(b, a, double(file_detrended_2D')).';

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
        file_data(:,:,:, idx2) = curr_file_data;
    end

    if bp == 1
        if detrended == 1
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(round(smoothing_kernel)), '_bp_d_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(round(smoothing_kernel)), '_bp_d_', file_name, '.gz']; % update cleaned_file
        else
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(round(smoothing_kernel)), '_bp_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(round(smoothing_kernel)), '_bp_', file_name, '.gz']; % update cleaned_file
        end
    else
        if detrended == 1
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(round(smoothing_kernel)), '_d_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(round(smoothing_kernel)), '_d_', file_name, '.gz']; % update cleaned_file
        else
            % write to data dir and relabel cleaned file:
            niftiwrite(cast(file_data, 'single'), [output_dir, '/s', num2str(round(smoothing_kernel)), '_', file_name], file_orig_data_info, "Compressed", true)
            cleaned_file = [output_dir, '/s', num2str(round(smoothing_kernel)), '_', file_name, '.gz']; % update cleaned_file
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