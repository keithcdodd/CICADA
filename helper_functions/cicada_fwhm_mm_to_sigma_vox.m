function [sigma_vox, sigma_mm, voxel_sizes_mm] = cicada_fwhm_mm_to_sigma_vox(fwhm_mm, voxel_sizes_mm)
%CICADA_FWHM_MM_TO_SIGMA_VOX Convert a physical Gaussian FWHM to voxel sigma.
%   [SIGMA_VOX, SIGMA_MM, VOXEL_SIZES_MM] =
%   CICADA_FWHM_MM_TO_SIGMA_VOX(FWHM_MM, VOXEL_SIZES_MM) validates a
%   nonnegative scalar FWHM and three positive spatial voxel dimensions.
%   SIGMA_VOX is a 1-by-3 vector suitable for IMGAUSSFILT3, preserving the
%   requested physical smoothing scale for anisotropic voxel geometry.

validateattributes(fwhm_mm, {'numeric'}, ...
    {'real','scalar','finite','nonnegative'}, mfilename, 'fwhm_mm');
validateattributes(voxel_sizes_mm, {'numeric'}, ...
    {'real','vector','numel',3,'finite','positive'}, ...
    mfilename, 'voxel_sizes_mm');

fwhm_mm = double(fwhm_mm);
voxel_sizes_mm = reshape(double(voxel_sizes_mm), 1, 3);
sigma_mm = fwhm_mm / (2 * sqrt(2 * log(2)));
sigma_vox = sigma_mm ./ voxel_sizes_mm;
end
