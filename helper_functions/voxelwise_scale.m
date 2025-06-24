function voxelwise_scaled_nii = voxelwise_scale(input_nii)
% voxel-wise scale scales each voxel time series to mean 100, this helps
% voxels be more comparable in group data.
%   output_nii = voxelwise_mean100_normalize(input_nii)
%
%   input_nii: Full path to input 4D .nii.gz file
%   output_nii: Full path to output scaled .nii.gz file (auto-named)

    % Ensure FSL is available
    if isempty(which('call_fsl'))
        error('call_fsl not found. Make sure FSL and call_fsl are configured.');
    end

    % Output path (same directory, suffix _mean100.nii.gz)
    [in_path, in_name, ~] = fileparts(input_nii);
    if endsWith(in_name, '.nii')  % handle .nii.gz double extension
        in_name = extractBefore(in_name, '.nii');
    end

    voxelwise_scaled_nii = fullfile(in_path, ['n_', in_name '.nii.gz']); % n for normalization

    % Temp files
    mean_img = fullfile(tempdir, [in_name '_mean.nii.gz']);
    norm_img = fullfile(tempdir, [in_name '_norm.nii.gz']);

    % Step 1: Compute voxel-wise mean over time
    call_fsl(['fslmaths ' input_nii ' -Tmean ' mean_img]);

    % Step 2: Divide input by mean (voxel-wise)
    call_fsl(['fslmaths ' input_nii ' -div ' mean_img ' ' norm_img]);

    % Step 3: Multiply by 100
    call_fsl(['fslmaths ' norm_img ' -mul 100 ' voxelwise_scaled_nii]);

    % Optional: clean up temp files
    delete(mean_img);
    delete(norm_img);
    delete(input_nii)
end
