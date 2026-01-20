function [refinedFuncMaskFile, refinedAnatMaskFile] = refine_masks(anatMaskFile, funcMaskFile, func4DFile)
% refine_masks
% Inputs:
%   anatMaskFile = path to anatomical mask (.nii.gz)
%   funcMaskFile = path to functional mask (.nii.gz) [not directly used, but input kept]
%   func4DFile   = path to functional 4D data (.nii.gz)
%
% Outputs:
%   refinedFuncMaskFile = path to refined functional mask (.nii.gz)
%   refinedAnatMaskFile = path to refined anatomical mask (.nii.gz)

%% ---- Load functional data (reference space) ----
funcInfo = niftiinfo(func4DFile);
funcData = double(niftiread(funcInfo));
meanFunc = mean(funcData, 4);

% Threshold at bottom 5% of nonzero voxels
nzVals = meanFunc(meanFunc > 0);
%thr = quantile(nzVals, 0.05);
thr = 0;
funcMask = meanFunc > thr;

%% ---- Load anatomical mask ----
anatInfo = niftiinfo(anatMaskFile);
anatMask = logical(niftiread(anatInfo));

%% ---- Resample anatomical mask into functional space ----
anatResampled = maybe_resample_mask_like(funcData(:,:,:,1), anatMask);

%% ---- Morphological operations ----
voxelSize = funcInfo.PixelDimensions; % use functional voxel size
se = createSphericalSE(3, voxelSize);

anatDilated3 = imdilate(imdilate(imdilate(anatResampled, se), se), se); % thrice dilated
refinedAnatMask = imerode(imerode(imerode(anatDilated3, se), se), se); % erode back to anatomical, just smoother

anatDilated2 = imdilate(imdilate(refinedAnatMask, se), se);

% Refined functional mask
refinedFuncMask = anatDilated2 & funcMask;

%% ---- Save results ----
[outPath, ~, ~] = fileparts(func4DFile);
refinedFuncMaskFile = fullfile(outPath, 'refined_funcmask.nii.gz');
refinedAnatMaskFile = fullfile(outPath, 'refined_anatmask.nii.gz');

funcOutInfo = funcInfo;
funcOutInfo.Datatype = 'uint8';
funcOutInfo.PixelDimensions = funcOutInfo.PixelDimensions(1:3);
funcOutInfo.ImageSize = funcOutInfo.ImageSize(1:3);
niftiwrite(uint8(refinedFuncMask), fullfile(outPath, 'refined_funcmask'), funcOutInfo, 'Compressed', true);

anatOutInfo = funcInfo; % save anatomical mask in functional space
anatOutInfo.Datatype = 'uint8';
anatOutInfo.PixelDimensions = funcOutInfo.PixelDimensions(1:3);
anatOutInfo.ImageSize = funcOutInfo.ImageSize(1:3);
niftiwrite(uint8(refinedAnatMask), fullfile(outPath, 'refined_anatmask'), anatOutInfo, 'Compressed', true);

end


function se = createSphericalSE(radius_mm, voxelSize)
% Create 3D spherical structuring element with physical radius in mm
radius_vox = radius_mm ./ voxelSize; 
[x, y, z] = ndgrid(-ceil(radius_vox(1)):ceil(radius_vox(1)), ...
                   -ceil(radius_vox(2)):ceil(radius_vox(2)), ...
                   -ceil(radius_vox(3)):ceil(radius_vox(3)));
dist = (x .* voxelSize(1)).^2 + (y .* voxelSize(2)).^2 + (z .* voxelSize(3)).^2;
se = strel(dist <= radius_mm^2);
end


function R = nifti2ref(info)
% Convert niftiinfo into imref3d object
dims = info.ImageSize;
pixdim = info.PixelDimensions;
R = imref3d(dims, ...
    pixdim(2) * (0.5:dims(2)-0.5), ...
    pixdim(1) * (0.5:dims(1)-0.5), ...
    pixdim(3) * (0.5:dims(3)-0.5));
end

function Vout = maybe_resample_mask_like(targetVol, sourceVol)
    % If size already matches, return as-is
    if isequal(size(targetVol,1), size(sourceVol,1)) && ...
       isequal(size(targetVol,2), size(sourceVol,2)) && ...
       isequal(size(targetVol,3), size(sourceVol,3))
        Vout = sourceVol;
        return;
    end

    % Resample with nearest-neighbor (preserve binary mask nature)
    try
        Vout = imresize3(sourceVol, size(targetVol), 'nearest');
    catch
        % Fallback manual nearest-neighbor
        [X,Y,Z] = ndgrid( linspace(1,size(sourceVol,1),size(targetVol,1)), ...
                          linspace(1,size(sourceVol,2),size(targetVol,2)), ...
                          linspace(1,size(sourceVol,3),size(targetVol,3)) );
        Vout = interpn(sourceVol, X, Y, Z, 'nearest', 0);
    end

    % Ensure strict binary output
    Vout = Vout > 0.5;
end
