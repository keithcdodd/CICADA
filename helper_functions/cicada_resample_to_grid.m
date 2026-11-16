function out = cicada_resample_to_grid(vol, srcInfo, refInfo, method)
% Affine-aware 3D resampling between grids in the SAME
% anatomical/world coordinate space.
%
% method: 'linear' for probabilities, 'nearest' for masks.

srcT = double(srcInfo.Transform.T);
refT = double(refInfo.Transform.T);
refSize = refInfo.ImageSize(1:3);

assert(all(isfinite(srcT(:))) && all(isfinite(refT(:))));
assert(isequal(size(srcT), [4 4]) && ...
       isequal(size(refT), [4 4]), ...
    'Source and reference affines must be 4x4.');

assert(rcond(srcT) > 1e-12, ...
    'Invalid source affine.');

assert(rcond(refT) > 1e-12, ...
    'Invalid reference affine.');

assert(isequal([size(vol,1), size(vol,2), size(vol,3)], ...
               srcInfo.ImageSize(1:3)), ...
    'Source volume dimensions do not match its NIfTI metadata.');

if isfield(srcInfo, 'SpaceUnits') && ...
        isfield(refInfo, 'SpaceUnits')

    assert(strcmpi(string(srcInfo.SpaceUnits), ...
                   string(refInfo.SpaceUnits)), ...
        'Source and reference spatial units differ.');
end

% MATLAB affine3d uses row-vector coordinates.
% Reference voxel -> world -> source voxel.
M = refT / srcT;

vol = double(vol);
out = zeros(refSize, 'double');

[i, j] = ndgrid(1:refSize(1), 1:refSize(2));

for k = 1:refSize(3)

    pts = [i(:), j(:), ...
           k * ones(numel(i), 1), ...
           ones(numel(i), 1)] * M;

    vals = interpn(vol, ...
        pts(:,1), pts(:,2), pts(:,3), method, 0);

    out(:,:,k) = reshape(vals, refSize(1:2));

end


end