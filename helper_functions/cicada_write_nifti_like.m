function output_path = cicada_write_nifti_like(data, reference_path, output_path, datatype)
% CICADA_WRITE_NIFTI_LIKE Write a 3-D NIfTI using reference geometry.
if nargin < 4 || isempty(datatype)
    datatype = class(data);
end
reference_path = char(reference_path);
output_path = char(output_path);
if ~isfile(reference_path)
    error('CICADA:MissingNIfTIReference', 'Reference NIfTI does not exist: %s', reference_path);
end
out_parent = fileparts(output_path);
if ~isempty(out_parent) && ~isfolder(out_parent)
    mkdir(out_parent);
end

out_data = cast(data, datatype);
if ndims(out_data) ~= 3
    error('CICADA:NIfTIOutputDimensions', 'cicada_write_nifti_like expects a 3-D array.');
end
info = niftiinfo(reference_path);
info.ImageSize = size(out_data);
if numel(info.PixelDimensions) > 3
    info.PixelDimensions = info.PixelDimensions(1:3);
end
info.Datatype = class(out_data);
if isfield(info, 'BitsPerPixel')
    info.BitsPerPixel = local_bits_per_pixel(class(out_data));
end
if isfield(info, 'MultiplicativeScaling')
    info.MultiplicativeScaling = 1;
end
if isfield(info, 'AdditiveOffset')
    info.AdditiveOffset = 0;
end

if endsWith(output_path, '.nii.gz')
    base_path = output_path(1:end-7);
    niftiwrite(out_data, base_path, info, 'Compressed', true);
elseif endsWith(output_path, '.nii')
    base_path = output_path(1:end-4);
    niftiwrite(out_data, base_path, info, 'Compressed', false);
else
    niftiwrite(out_data, output_path, info, 'Compressed', false);
end
if ~isfile(output_path)
    error('CICADA:NIfTIWriteFailure', 'Expected NIfTI was not created: %s', output_path);
end
end

function bits = local_bits_per_pixel(datatype)
switch datatype
    case {'uint8','int8','logical'}
        bits = 8;
    case {'uint16','int16'}
        bits = 16;
    case {'uint32','int32','single'}
        bits = 32;
    case {'uint64','int64','double'}
        bits = 64;
    otherwise
        error('CICADA:UnsupportedNIfTIDatatype', ...
            'Unsupported output datatype: %s', datatype);
end
end
