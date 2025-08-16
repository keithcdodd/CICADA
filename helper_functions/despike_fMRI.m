function [despiked_file, spikeIdxPerVoxel] = despike_fmri(funcFile, gmProbFile, robustZ_thresh)
% robustZ_thresh is typically 3-5. 3 is usually good, 5 is more
% conservative.

hp_cutoff = 0.008;
%% Resample GM probability map to functional resolution using FSL
[funcPath, name, ext] = fileparts(funcFile);
if strcmpi(ext, '.gz')
    % strip the .gz and check again
    [~, name, ext2] = fileparts(name);
    name = [name]; % name without extension
end

resampledGM = fullfile(gmProbFile);

% Load GM and functional info
funcInfo = niftiinfo(funcFile);
rawFunc = double(niftiread(funcInfo));
Y4D = rawFunc * funcInfo.MultiplicativeScaling + funcInfo.AdditiveOffset;
szFunc = size(Y4D);

gmInfo = niftiinfo(resampledGM);
rawGM = double(niftiread(gmInfo));
gmProb = rawGM * gmInfo.MultiplicativeScaling + gmInfo.AdditiveOffset;

% Define reference objects for spatial mapping
Rfunc = imref3d(szFunc(1:3), ...
    funcInfo.PixelDimensions(1), funcInfo.PixelDimensions(2), funcInfo.PixelDimensions(3));
Rgm = imref3d(size(gmProb), ...
    gmInfo.PixelDimensions(1), gmInfo.PixelDimensions(2), gmInfo.PixelDimensions(3));

% Compute grid in GM space corresponding to functional voxels
[xF, yF, zF] = ndgrid( (0:szFunc(1)-1)*Rfunc.PixelExtentInWorldX, ...
                        (0:szFunc(2)-1)*Rfunc.PixelExtentInWorldY, ...
                        (0:szFunc(3)-1)*Rfunc.PixelExtentInWorldZ );

% Map functional world coordinates to GM voxel coordinates
xGM = xF / Rgm.PixelExtentInWorldX + 1;
yGM = yF / Rgm.PixelExtentInWorldY + 1;
zGM = zF / Rgm.PixelExtentInWorldZ + 1;

% Trilinear interpolation
resampledGM = interp3(double(gmProb), yGM, xGM, zGM, 'linear', 0); 
% Note interp3 uses (Y,X,Z) indexing

% Threshold to create mask
gmMask = resampledGM > 0.5;


% Get TR from NIfTI header
TR = funcInfo.PixelDimensions(4);
if TR <= 0
    error('TR could not be determined from NIfTI header.');
end

%% Parameters
sz = size(Y4D);
T = sz(4);
Fs = 1/TR;

% High-pass filter design
Wn = hp_cutoff / (Fs/2);
[b,a] = butter(2, Wn, 'high');

% Output init
Y_out = Y4D;
gmIdx = find(gmMask);
nVox = numel(gmIdx);
spikeIdxPerVoxel = cell(nVox,1);
mad1 = @(x) median(abs(x - median(x))) + eps;

altered_voxels_before = [];
altered_voxels_after  = [];

%% Loop over GM voxels
for k = 1:nVox
    [ix,iy,iz] = ind2sub(sz(1:3), gmIdx(k));
    ts = squeeze(Y4D(ix,iy,iz,:));

    if any(ts == 0)
        continue;
    end

    ts_hp = filtfilt(b,a,ts);
    r = ts_hp - median(ts_hp);
    rz = abs(r) / (1.4826 * mad1(r));

    spikeIdx = find(rz > robustZ_thresh);
    if isempty(spikeIdx)
        continue;
    end

    repl = ts;
    bad  = false(T,1);
    bad(spikeIdx) = true;

    tmp = repl;
    tmp(bad) = NaN;
    tmp = fillmissing(tmp, 'linear');

    if isnan(tmp(1)), tmp(1) = tmp(find(~isnan(tmp),1,'first')); end
    if isnan(tmp(end)), tmp(end) = tmp(find(~isnan(tmp),1,'last')); end

    repl(bad) = tmp(bad);
    Y_out(ix,iy,iz,:) = repl;

    spikeIdxPerVoxel{k} = spikeIdx(:);

    altered_voxels_before = [altered_voxels_before, filtfilt(b,a,ts)]; %#ok<AGROW>
    altered_voxels_after  = [altered_voxels_after,  filtfilt(b,a,repl)]; %#ok<AGROW>
end

%% QC plot (save only)
if ~isempty(altered_voxels_before)
    mean_before = mean(altered_voxels_before, 2, 'omitnan');
    mean_after  = mean(altered_voxels_after,  2, 'omitnan');

    % Create invisible figure
    hFig = figure('Visible','off');
    plot((0:T-1)*TR, mean_before, 'r', 'LineWidth', 1.5); hold on;
    plot((0:T-1)*TR, mean_after, 'b', 'LineWidth', 1.5);
    xlabel('Time (s)');
    ylabel('Signal (a.u., HPF)');
    legend({'Before despike','After despike'}, 'Location', 'best');
    title('Mean HPF time series of altered voxels');
    grid on;

    % Save figure as PNG
    saveas(hFig, fullfile(fileparts(funcFile), [name, '_despike_QC.png']));

    % Close figure to free memory
    close(hFig);
end


%% Save output NIfTI (rescale to match original storage)
scaledBack = (Y_out - funcInfo.AdditiveOffset) ./ funcInfo.MultiplicativeScaling;

% Match header datatype to stored data
%funcInfo.Datatype = class(rawFunc);
%funcInfo.BitsPerPixel = 8 * numel(typecast(cast(0,class(rawFunc)),'uint8'));

outNii = fullfile(fileparts(funcFile), [name, '_despiked']);
niftiwrite(cast(scaledBack, 'int16'), outNii, funcInfo, "Compressed", true);

despiked_file = [outNii, '.nii.gz'];
end
