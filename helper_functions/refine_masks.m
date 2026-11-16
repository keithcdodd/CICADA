function refinedFuncMaskFile = refine_masks(anatMaskFile, funcMaskFile, func4DFile)
% Reasonable to include for most data though is often likely not needed
% The idea is that sometimes fMRIPrep or other functions are too aggressive with functional masks
% This then refines the original funcmask to be less aggressive
% refine_masks
% Inputs:
%   funcMaskFile = path to functional mask (.nii.gz)
%                  Used to constrain valid functional coverage.
%   func4DFile   = path to functional 4D data (.nii.gz)
% Outputs:
%   refinedFuncMaskFile = path to refined functional mask (.nii.gz)

%% ---- Load functional data (reference space) ----
funcInfo = niftiinfo(func4DFile);

% Compute mean BOLD from stored voxel values.
funcData = single(niftiread(funcInfo));
meanFunc = double(mean(funcData, 4));
clear funcData

% Apply NIfTI intensity scaling.
slope = double(funcInfo.MultiplicativeScaling);
offset = double(funcInfo.AdditiveOffset);

assert(isfinite(slope), 'Invalid NIfTI scaling slope.');

if slope ~= 0
    assert(isfinite(offset), 'Invalid NIfTI scaling intercept.');
    meanFunc = meanFunc .* slope + offset;
end

funcMaskInfo = niftiinfo(funcMaskFile);
inputFuncMask = niftiread(funcMaskInfo) > 0;

assert(isequal(funcMaskInfo.ImageSize(1:3), ...
               funcInfo.ImageSize(1:3)), ...
    'Functional mask and BOLD dimensions differ.');

assert(numel(funcMaskInfo.ImageSize) == 3, ...
    'Functional mask must be a 3D image.');

assert(max(abs(funcMaskInfo.Transform.T(:) - ...
               funcInfo.Transform.T(:))) < 1e-4, ...
    'Functional mask and BOLD affines differ.');

% Identify voxels containing finite, nonzero-mean BOLD data.
% Do not require membership in the original fMRIPrep mask.
validBOLD = isfinite(meanFunc) & (meanFunc > 0);

%% ---- Load anatomical mask ----
anatInfo = niftiinfo(anatMaskFile);
anatMask = logical(niftiread(anatInfo));

%% ---- Resample anatomical mask into functional space ----
anatResampled = cicada_resample_to_grid( ...
    double(anatMask), anatInfo, funcInfo, 'nearest') > 0.5;

%% ---- Repair missing functional brain coverage ----

% Anatomical brain mask resampled into functional space.
refinedAnatMask = anatResampled;

% Identify brain voxels excluded from the original functional mask.
recoveryCandidates = refinedAnatMask & ~inputFuncMask;

%% ---- FSL robust intensity criterion for recovered voxels ----

% Create anatomical-masked mean BOLD using physical intensities.
maskedMean = zeros(size(meanFunc), 'single');

validAnat = refinedAnatMask & ...
            isfinite(meanFunc) & (meanFunc > 0);

assert(any(validAnat(:)), ...
    'No positive BOLD intensities within the anatomical mask.');

maskedMean(validAnat) = single(meanFunc(validAnat));

% Create a temporary directory for FSL calculations.
tmpDir = tempname;
mkdir(tmpDir);
cleanupTmp = onCleanup(@() rmdir(tmpDir, 's'));

% Save the mean image with 3D functional-space geometry.
tmpInfo = funcInfo;
tmpInfo.ImageSize = funcInfo.ImageSize(1:3);
tmpInfo.PixelDimensions = funcInfo.PixelDimensions(1:3);
tmpInfo.Datatype = 'single';
tmpInfo.BitsPerPixel = 32;
tmpInfo.MultiplicativeScaling = 1;
tmpInfo.AdditiveOffset = 0;

meanFile = fullfile(tmpDir, 'mean_anatmasked.nii.gz');
robustFile = fullfile(tmpDir, 'robust5_mask.nii.gz');

niftiwrite(maskedMean, ...
    fullfile(tmpDir, 'mean_anatmasked'), ...
    tmpInfo, 'Compressed', true);

% Apply FSL's robust 5% threshold to positive mean EPI intensities.
cmd = sprintf('fslmaths "%s" -thrP 5 -bin "%s"', ...
    meanFile, robustFile);

[status, output] = system(cmd);

if status ~= 0
    error('FSL robust thresholding failed: %s', output);
end

robustSignal = niftiread(robustFile) > 0;

assert(isequal(size(robustSignal), size(refinedAnatMask)), ...
    'FSL robust mask dimensions do not match functional space.');

% Recover only anatomically supported voxels passing the threshold.
recoveredVoxels = recoveryCandidates & robustSignal;

% Preserve the original mask and add supported missing voxels.
refinedFuncMask = inputFuncMask | recoveredVoxels;

% Quantify reconstruction.
nOriginal  = nnz(inputFuncMask);
nCandidate = nnz(recoveryCandidates);
nAdded     = nnz(recoveredVoxels);
nFinal     = nnz(refinedFuncMask);

fprintf('\nFunctional mask reconstruction:\n');

fprintf('  Original mask:       %d voxels\n', nOriginal);
fprintf('  Anatomical gaps:     %d voxels\n', nCandidate);
fprintf('  Recovered voxels:    %d voxels\n', nAdded);
fprintf('  Final mask:          %d voxels\n', nFinal);

fprintf('  Increase:            %.2f%%\n', ...
    100 * nAdded / max(nOriginal, 1));

% Compare recovered voxel intensities with existing brain coverage.
referenceVoxels = inputFuncMask & refinedAnatMask & validBOLD;

if nAdded > 0 && any(referenceVoxels(:))

    referenceMedian = median(meanFunc(referenceVoxels));
    recoveredMedian = median(meanFunc(recoveredVoxels));

    fprintf('  Median recovered/reference EPI intensity: %.3f\n', ...
        recoveredMedian / referenceMedian);
end

%% ---- Save results ----
[outPath, ~, ~] = fileparts(func4DFile);
refinedFuncMaskFile = fullfile(outPath, 'refined_funcmask.nii.gz');

funcOutInfo = funcInfo;
funcOutInfo.Datatype = 'uint8';
funcOutInfo.BitsPerPixel = 8;
funcOutInfo.MultiplicativeScaling = 1;
funcOutInfo.AdditiveOffset = 0;
funcOutInfo.PixelDimensions = funcOutInfo.PixelDimensions(1:3);
funcOutInfo.ImageSize = funcOutInfo.ImageSize(1:3);
niftiwrite(uint8(refinedFuncMask), fullfile(outPath, 'refined_funcmask'), funcOutInfo, 'Compressed', true);

% anatOutInfo = funcInfo; % save anatomical mask in functional space
% anatOutInfo.Datatype = 'uint8';
% anatOutInfo.BitsPerPixel = 8;
% anatOutInfo.MultiplicativeScaling = 1;
% anatOutInfo.AdditiveOffset = 0;
% anatOutInfo.PixelDimensions = funcOutInfo.PixelDimensions(1:3);
% anatOutInfo.ImageSize = funcOutInfo.ImageSize(1:3);
% niftiwrite(uint8(refinedAnatMask), fullfile(outPath, 'refined_anatmask'), anatOutInfo, 'Compressed', true);

end