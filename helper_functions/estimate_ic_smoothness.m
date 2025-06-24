function fwhm_avg = estimate_ic_smoothness(ic_1d, func_mask_3d_nii, voxel_dim_3D, voxel_size_3D)
% Estimate FWHM smoothness for a single subject-specific IC map using FSL
%
% Inputs:
%   - ic_1d: [V × 1] vector of a single dual-regression IC map (V = voxels)
%   - func_mask_3d_nii: 3D binary mask .nii.gz file
%   - voxel_dim_3D: [nx, ny, nz] dimensions of the 3D volume
%   - voxel_size_3D: [vx, vy, vz] voxel sizes in mm
%
% Output:
%   - fwhm_avg: scalar average FWHM across x, y, z (in mm)

    [nx, ny, nz] = deal(voxel_dim_3D(1), voxel_dim_3D(2), voxel_dim_3D(3));
    nVox = nx * ny * nz;
    
    if length(ic_1d) ~= nVox
        error('IC vector length does not match volume size');
    end

    % Reshape IC to 3D
    ic_3d = reshape(ic_1d, [nx, ny, nz]);

    % Prepare NIfTI info
    funcmask_info = niftiinfo(func_mask_3d_nii);
    funcmask_info.PixelDimensions = voxel_size_3D;
    funcmask_info.ImageSize = [nx, ny, nz];
    funcmask_info.DataType = 'single';

    % Temporary filenames
    temp_prefix = tempname;
    ic_file = temp_prefix + "_ic.nii.gz";
    mask_file = func_mask_3d_nii;

    % Save IC map
    niftiwrite(single(ic_3d), [temp_prefix, '_ic'], funcmask_info, 'Compressed', true);

    % Call FSL smoothest
    fsl_cmd = sprintf('smoothest -z %s -m %s', ic_file, mask_file);
    [status, output] = call_fsl(fsl_cmd);

    % Parse output
    if status == 0
        % Extract FWHMmm values using regexp
        tokens = regexp(output, 'FWHMmm ([\d\.]+) ([\d\.]+) ([\d\.]+)', 'tokens');
        
        % Convert from cell to numeric
        if ~isempty(tokens)
            fwhm_vals = str2double(tokens{1});
            fwhm_avg = mean(fwhm_vals);
        else
            warning('FWHMmm values not found in the output');
            fwhm_avg = NaN;
        end
    else
        warning('FSL smoothest failed:\n%s', output);
        fwhm_avg = NaN;
    end

    % Clean up
    if exist(ic_file, 'file'), delete(ic_file); end
end
