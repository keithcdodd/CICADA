function cicada_assert_nifti_geometry(reference_path, varargin)
% CICADA_ASSERT_NIFTI_GEOMETRY Require matching 3-D size, voxel size, and affine.
reference_path = char(reference_path);
if ~isfile(reference_path)
    error('CICADA:MissingGeometryReference', 'Reference NIfTI does not exist: %s', reference_path);
end
ref = niftiinfo(reference_path);
ref_size = double(ref.ImageSize(1:min(3,numel(ref.ImageSize))));
ref_pix = double(ref.PixelDimensions(1:min(3,numel(ref.PixelDimensions))));
ref_affine = local_affine(ref);
for i = 1:numel(varargin)
    candidate = char(varargin{i});
    if ~isfile(candidate)
        error('CICADA:MissingGeometryCandidate', 'NIfTI does not exist: %s', candidate);
    end
    info = niftiinfo(candidate);
    curr_size = double(info.ImageSize(1:min(3,numel(info.ImageSize))));
    curr_pix = double(info.PixelDimensions(1:min(3,numel(info.PixelDimensions))));
    curr_affine = local_affine(info);
    if ~isequal(ref_size, curr_size) || numel(ref_pix) ~= numel(curr_pix) || ...
            any(abs(ref_pix-curr_pix) > 1e-5) || ...
            ~isequal(size(ref_affine), size(curr_affine)) || ...
            any(abs(ref_affine(:)-curr_affine(:)) > 1e-4)
        error('CICADA:GeometryMismatch', ...
            'NIfTI geometry differs between %s and %s.', reference_path, candidate);
    end
end
end

function affine = local_affine(info)
affine = eye(4);
if ~isfield(info, 'Transform') || isempty(info.Transform)
    return
end
transform = info.Transform;
if isstruct(transform) && isfield(transform, 'T')
    affine = double(transform.T);
elseif isobject(transform) && isprop(transform, 'T')
    affine = double(transform.T);
end
end
