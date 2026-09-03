function [despiked_file, spikeIdxPerVoxel, qc] = despike_fMRI( ...
    funcFile, gmProbFile, cut2, outputDir)
% DESPIKE_FMRI_GM AFNI-inspired (afni.nimh.nih.gov) GM-only despiking for CICADA.
%
% This function follows the documented AFNI 3dDespike / -NEW + -localedit
% concepts, with one deliberate CICADA-specific spatial restriction: only GM
% voxels (GM probability > 0.5) can be edited.
%
% It is AFNI-inspired, not claimed numerically identical to AFNI.
%
% The idea: Very large spikes in data can make it difficult for MELODIC to
% separate signal ICs from noise ICs. Despiking the largest spikes only
% in GM may help better differentiate signal vs noise ICs (shown in small
% % sample size testing). Only GM voxels are affected so as to help limit
% unintended consequences with CICADA IC classification.
%
% Fit method:
%   - Before fitting, use AFNI -NEW's documented 9-point local-median
%     presquash concept: compute a 9-point local median and local MAD,
%     threshold at 6.789 * median(local MAD), and replace only those
%     preliminary extreme points by the local median FOR FITTING ONLY.
%   - Fit the quadratic + Fourier basis by least squares to that temporarily
%     presquashed copy. The original time series remains the data being
%     evaluated and edited.
%
% Final spike/edit method:
%   - sigma = sqrt(pi/2) * median(abs(original - fitted baseline))
%   - spikiness s = residual / sigma
%   - for |s| >= cut2, replace the original sample with the mean of the
%     nearest preceding and following non-spike samples
%   - outside GM: exactly unchanged
%
% CICADA default:
%   cut2 = 6
%
% This threshold was selected during CICADA development to restrict editing
% to extreme residual excursions while providing useful ICA preconditioning.
% Other positive cut2 values remain supported for explicit user testing.
%
% Despiking is voxelwise and does NOT censor or remove entire volumes.
% Volume-level edit burden is exported for QC and optional downstream review.

%% Input checks

if ~isfile(funcFile)
    error('Functional file not found: %s', funcFile);
end

if ~isfile(gmProbFile)
    error('GM probability file not found: %s', gmProbFile);
end

if ~exist('cut2', 'var') || isempty(cut2)
    cut2 = 6;
end

if ~isscalar(cut2) || ~isfinite(cut2) || cut2 <= 0
    error('cut2 must be a finite positive scalar.');
end

%% Output basename

[funcPath, name, ext] = fileparts(funcFile);

if strcmpi(ext, '.gz')
    [~, name, ext2] = fileparts(name);
    if ~strcmpi(ext2, '.nii')
        error('Expected a .nii or .nii.gz functional file.');
    end
elseif ~strcmpi(ext, '.nii')
    error('Expected a .nii or .nii.gz functional file.');
end

if ~exist('outputDir', 'var') || isempty(outputDir)
    outputDir = funcPath;
end

if ~ischar(outputDir) && ~isstring(outputDir)
    error('outputDir must be a character array or string.');
end

outputDir = char(outputDir);

if ~isfolder(outputDir)
    mkdir(outputDir);
end

%% Load functional data while preserving storage information

funcInfo = niftiinfo(funcFile);
storedFunc = niftiread(funcInfo);
storedClass = class(storedFunc);

rawFunc = double(storedFunc);

scale = funcInfo.MultiplicativeScaling;
offset = funcInfo.AdditiveOffset;

if isempty(scale) || ~isfinite(scale) || scale == 0
    scale = 1;
end

if isempty(offset) || ~isfinite(offset)
    offset = 0;
end

Y4D = rawFunc .* scale + offset;
szFunc = size(Y4D);

if numel(szFunc) < 4 || szFunc(4) < 15
    error('Functional input must be 4-D with at least 15 time points.');
end

T = szFunc(4);

%% Load GM probability map

gmInfo = niftiinfo(gmProbFile);
storedGM = niftiread(gmInfo);
rawGM = double(storedGM);

gmScale = gmInfo.MultiplicativeScaling;
gmOffset = gmInfo.AdditiveOffset;

if isempty(gmScale) || ~isfinite(gmScale) || gmScale == 0
    gmScale = 1;
end

if isempty(gmOffset) || ~isfinite(gmOffset)
    gmOffset = 0;
end

gmProb = squeeze(rawGM .* gmScale + gmOffset);

if ndims(gmProb) ~= 3
    error('GM probability input must be a 3-D volume.');
end

%% Affine-aware resampling of GM probability map to functional grid

if ~isfield(funcInfo, 'Transform') || isempty(funcInfo.Transform) || ...
        ~isfield(gmInfo, 'Transform') || isempty(gmInfo.Transform)
    error('Functional and GM NIfTI files must contain valid spatial transforms.');
end

sameSpatialSize = isequal(size(gmProb), szFunc(1:3));

testPoints = [1 1 1; 2 1 1; 1 2 1; 1 1 2];

try
    funcWorld = transformPointsForward(funcInfo.Transform, testPoints);
    gmWorld = transformPointsForward(gmInfo.Transform, testPoints);
    transformDifference = max(abs(funcWorld(:) - gmWorld(:)));
catch
    transformDifference = Inf;
end

if sameSpatialSize && transformDifference < 1e-6
    resampledGM = gmProb;
else
    nRows   = szFunc(1);
    nCols   = szFunc(2);
    nSlices = szFunc(3);

    [rowF, colF, sliceF] = ndgrid(1:nRows, 1:nCols, 1:nSlices);

    [xWorld, yWorld, zWorld] = transformPointsForward( ...
        funcInfo.Transform, colF, rowF, sliceF);

    [colGM, rowGM, sliceGM] = transformPointsInverse( ...
        gmInfo.Transform, xWorld, yWorld, zWorld);

    resampledGM = interp3( ...
        double(gmProb), colGM, rowGM, sliceGM, 'linear', 0);
end

gmMask = resampledGM > 0.5;

if ~any(gmMask(:))
    error(['GM mask contains no voxels after alignment/resampling. ', ...
        'Check the GM probability map and functional-image geometry.']);
end

%% TR, for QC plotting

if numel(funcInfo.PixelDimensions) < 4
    error('TR could not be determined from NIfTI header.');
end

TR = funcInfo.PixelDimensions(4);

if ~isfinite(TR) || TR <= 0
    error('TR could not be determined from NIfTI header.');
end

%% Build temporal basis

corder = round(T / 30);
nRef = 2 * corder + 3;

if nRef >= T
    error('Temporal basis too large for T=%d (corder=%d, nRef=%d).', ...
        T, corder, nRef);
end

t = (0:T-1)';
tm = 0.5 * (T - 1);
fac = 2 / T;
u = (t - tm) .* fac;

X = [ones(T,1), u, u.^2];

for kk = 1:corder
    fq = (2*pi*kk) / T;
    X = [X, sin(fq*t), cos(fq*t)]; %#ok<AGROW>
end

fprintf(['Despiking Parameters: T=%d, TR=%.6g s, ', ...
    'corder=%d, regressors=%d, cut2=%.6g\n'], ...
    T, TR, corder, size(X,2), cut2);

%% Initialize output and QC

Y_out = Y4D;

gmIdx = find(gmMask);
nVox = numel(gmIdx);

if nargout > 1
    spikeIdxPerVoxel = cell(nVox, 1);
else
    spikeIdxPerVoxel = [];
end

flagsPerVolume = zeros(T,1);

% Aggregate QC measures.
sumAbsResidualBefore = zeros(T,1);
sumAbsResidualAfter  = zeros(T,1);
sumAbsChange         = zeros(T,1);
sumSqChange          = zeros(T,1);
maxAbsChangePerVolume = zeros(T,1);

% Keep one concrete example: the voxel containing the largest actual
% localedit-induced signal change anywhere in GM.
exampleMaxAbsChange = -Inf;
exampleVoxelSub     = [];
exampleTs           = [];
exampleRepl         = [];
exampleSpikeIdx     = [];

nAlteredVoxels = 0;
nAlteredSamples = 0;
nZeroSkipped = 0;
nFitFailed = 0;
maxAbsSpikinessAltered = 0;

% Precompute 9-point window indices. Shift the first/last
% windows so every location uses a full 9-point window.
if T < 9
    error('Fitting requires at least 9 time points.');
end

windowStart = (1:T)' - 4;
windowStart = max(windowStart, 1);
windowStart = min(windowStart, T - 8);
afniWindowIdx = windowStart + (0:8);

%% Loop over GM voxels

fprintf('Processing %d GM voxels...\n', nVox);

for k = 1:nVox

    if mod(k,5000) == 0
        fprintf('  Processed %d / %d GM voxels...\n', k, nVox);
    end

    [ix, iy, iz] = ind2sub(szFunc(1:3), gmIdx(k));
    ts = squeeze(Y4D(ix, iy, iz, :));

    % CICADA-specific conservative guard retained from current function.
    if any(ts == 0) || any(~isfinite(ts))
        nZeroSkipped = nZeroSkipped + 1;
        continue;
    end

    try
        % Simplistic local-median presquash is used only
        % to make the smooth baseline fit insensitive to catastrophic spikes.
        fitInput = cicada_new_prefit(ts, afniWindowIdx);
        beta = X \ fitInput;
    catch
        nFitFailed = nFitFailed + 1;
        continue;
    end

    fitTs = X * beta;
    residual = ts - fitTs;

    % Scale definition.
    sigma = sqrt(pi/2) * median(abs(residual));

    if ~isfinite(sigma) || sigma <= 0
        continue;
    end

    spikiness = residual ./ sigma;
    isSpike = abs(spikiness) >= cut2;
    spikeIdx = find(isSpike);

    if isempty(spikeIdx)
        continue;
    end

    repl = ts;

    % Replacement local edit: Each flagged sample uses the nearest
    % non-spike sample on each side. Consecutive spike runs therefore share
    % the same bounding clean-pair mean. At temporal edges, use the nearest
    % available non-spike sample; only the pathological no-neighbor case falls
    % back to the fitted baseline.
    for jj = 1:numel(spikeIdx)

        iv = spikeIdx(jj);

        prevIdx = find(~isSpike(1:iv-1), 1, 'last');
        nextRelative = find(~isSpike(iv+1:end), 1, 'first');

        if isempty(nextRelative)
            nextIdx = [];
        else
            nextIdx = iv + nextRelative;
        end

        if ~isempty(prevIdx) && ~isempty(nextIdx)
            repl(iv) = 0.5 * (ts(prevIdx) + ts(nextIdx));
        elseif ~isempty(nextIdx)
            repl(iv) = ts(nextIdx);
        elseif ~isempty(prevIdx)
            repl(iv) = ts(prevIdx);
        else
            repl(iv) = fitTs(iv);
        end
    end

    Y_out(ix, iy, iz, :) = repl;

    if nargout > 1
        spikeIdxPerVoxel{k} = spikeIdx(:);
    end

    flagsPerVolume(spikeIdx) = flagsPerVolume(spikeIdx) + 1;

    % Actual change produced by localedit.
    delta = repl - ts;

    % Absolute residuals avoid cancellation across voxels.
    sumAbsResidualBefore = sumAbsResidualBefore + abs(residual);
    sumAbsResidualAfter  = sumAbsResidualAfter  + abs(repl - fitTs);

    % Magnitude of the actual modification to the BOLD data.
    sumAbsChange = sumAbsChange + abs(delta);
    sumSqChange  = sumSqChange  + delta.^2;
    maxAbsChangePerVolume = max(maxAbsChangePerVolume, abs(delta));

    % Save the voxel containing the largest single localedit change so the
    % QC figure can show an interpretable original-vs-despiked example.
    [localMaxAbsChange, ~] = max(abs(delta));

    if localMaxAbsChange > exampleMaxAbsChange
        exampleMaxAbsChange = localMaxAbsChange;
        exampleVoxelSub = [ix iy iz];
        exampleTs       = ts;
        exampleRepl     = repl;
        exampleSpikeIdx = spikeIdx(:);
    end

    nAlteredVoxels = nAlteredVoxels + 1;
    nAlteredSamples = nAlteredSamples + numel(spikeIdx);

    maxAbsSpikinessAltered = max( ...
        maxAbsSpikinessAltered, max(abs(spikiness(spikeIdx))));
end

if nAlteredVoxels == 0
    exampleMaxAbsChange = 0;
end

%% Summary

fprintf(['Summary: %d/%d GM voxels altered; ', ...
    '%d voxel-time samples replaced; cut2=%.6g.\n'], ...
    nAlteredVoxels, nVox, nAlteredSamples, cut2);

fprintf('  Time series skipped for zero/nonfinite values: %d\n', nZeroSkipped);
fprintf('  Baseline-fit failures: %d\n', nFitFailed);

%% QC structure

qc = struct();
qc.method = 'Local Edit, GM only';
qc.cut2 = cut2;
qc.gm_probability_threshold = 0.5;
qc.TR = TR;
qc.T = T;
qc.corder = corder;
qc.n_temporal_regressors = size(X,2);
qc.n_gm_voxels = nVox;
qc.n_altered_voxels = nAlteredVoxels;
qc.n_altered_samples = nAlteredSamples;
qc.fraction_gm_voxel_time_samples = nAlteredSamples / (nVox * T);
qc.flags_per_volume = flagsPerVolume;
qc.volumes_with_edits = nnz(flagsPerVolume);
qc.max_edits_one_volume = max(flagsPerVolume);
qc.edited_fraction_per_volume = flagsPerVolume ./ nVox;
qc.edited_volume_indices = find(flagsPerVolume > 0);
qc.n_zero_or_nonfinite_skipped = nZeroSkipped;
qc.fit_method = '9-point presquash + least squares';
qc.n_fit_failed = nFitFailed;
qc.max_abs_spikiness_among_edited_samples = maxAbsSpikinessAltered;
qc.max_abs_change_per_volume = maxAbsChangePerVolume;
qc.volume_censoring_applied = false;
qc.volume_qc_informational_only = true;
qc.motion_confounds_recommendation = [ ...
    'For CICADA IC classification, retain FD/DVARS from the original ', ...
    'preprocessing confounds file; do not recompute them from despiked data.'];
qc.replacement_rule = [ ...
    'Each flagged voxel-time sample is replaced by the mean of the ', ...
    'nearest preceding and following non-spike samples; at temporal ', ...
    'boundaries, the nearest available non-spike sample is used.'];
qc.boundary_rule = [ ...
    'The 9-point -NEW-inspired prefit window is shifted at temporal ', ...
    'edges so that all locations use a full 9-point window.'];

if nAlteredVoxels > 0
    qc.mean_abs_residual_before_altered_voxels = ...
        sumAbsResidualBefore ./ nAlteredVoxels;

    qc.mean_abs_residual_after_altered_voxels = ...
        sumAbsResidualAfter ./ nAlteredVoxels;
else
    qc.mean_abs_residual_before_altered_voxels = zeros(T,1);
    qc.mean_abs_residual_after_altered_voxels  = zeros(T,1);
end

% These use ALL GM voxels as the denominator. Voxels that were not
% edited therefore correctly contribute zero change.
qc.mean_abs_change_all_gm = sumAbsChange ./ nVox;
qc.rms_change_all_gm      = sqrt(sumSqChange ./ nVox);

qc.max_abs_signal_change = exampleMaxAbsChange;
qc.example_voxel_ijk     = exampleVoxelSub;
qc.example_spike_indices = exampleSpikeIdx;

%% Volume-level despiking QC table
%
% This table is informational only. CICADA does not automatically censor
% volumes based on despiking burden. Motion/DVARS confounds used elsewhere
% in CICADA should remain those from the original preprocessing pipeline.

cutTag = cicada_numeric_tag(cut2);

volumeOneBased = (1:T)';
volumeZeroBased = (0:T-1)';
timeSecQC = volumeZeroBased .* TR;

editedGMVoxels = flagsPerVolume;
editedGMFraction = flagsPerVolume ./ nVox;
meanAbsChangeAllGM = sumAbsChange ./ nVox;
rmsChangeAllGM = sqrt(sumSqChange ./ nVox);

volumeQC = table( ...
    volumeOneBased, ...
    volumeZeroBased, ...
    timeSecQC, ...
    editedGMVoxels, ...
    editedGMFraction, ...
    meanAbsChangeAllGM, ...
    rmsChangeAllGM, ...
    maxAbsChangePerVolume, ...
    'VariableNames', { ...
        'VolumeOneBased', ...
        'VolumeZeroBased', ...
        'TimeSeconds', ...
        'EditedGMVoxels', ...
        'EditedGMFraction', ...
        'MeanAbsChangeAllGM', ...
        'RMSChangeAllGM', ...
        'MaxAbsChangeGM'});

volumeQcFile = fullfile( ...
    outputDir, ...
    [name, '_localGM_c2-', cutTag, '_volumeQC.tsv']);

writetable(volumeQC, volumeQcFile, ...
    'FileType','text', ...
    'Delimiter','\t');

qc.volume_qc_file = volumeQcFile;

fprintf('Volume-level despiking QC: %s\n', volumeQcFile);

%% Enhanced QC plot

if nAlteredVoxels > 0

    timeSec = (0:T-1)' .* TR;

    meanAbsResidualBefore = ...
        sumAbsResidualBefore ./ nAlteredVoxels;

    meanAbsResidualAfter = ...
        sumAbsResidualAfter ./ nAlteredVoxels;

    hFig = figure( ...
        'Visible','off', ...
        'Color','w', ...
        'Position',[100 100 1200 950]);

    tl = tiledlayout(hFig,4,1, ...
        'TileSpacing','compact', ...
        'Padding','compact');


    % ------------------------------------------------------------
    % 1. How many GM voxels were edited at each volume?
    % ------------------------------------------------------------
    nexttile

    bar(timeSec, flagsPerVolume, 1);

    ylabel('Edited GM voxels');
    title('Localedit burden by volume');
    xlim([timeSec(1) timeSec(end)]);
    grid on


    % ------------------------------------------------------------
    % 2. Absolute residual before vs after localedit
    % ------------------------------------------------------------
    nexttile

    plot(timeSec, meanAbsResidualBefore, ...
        'LineWidth',1.25);

    hold on

    plot(timeSec, meanAbsResidualAfter, ...
        'LineWidth',1.25);

    ylabel('Mean |residual| (a.u.)');

    legend( ...
        {'Before localedit','After localedit'}, ...
        'Location','best');

    title('Mean absolute residual across altered GM voxels');
    xlim([timeSec(1) timeSec(end)]);
    grid on


    % ------------------------------------------------------------
    % 3. Magnitude of the actual data modification
    % ------------------------------------------------------------
    nexttile

    plot(timeSec, rmsChangeAllGM, ...
        'LineWidth',1.25);

    ylabel('RMS change (a.u.)');
    title('RMS localedit-induced signal change across all GM voxels');
    xlim([timeSec(1) timeSec(end)]);
    grid on


    % ------------------------------------------------------------
    % 4. Concrete example: voxel with largest actual edit
    % ------------------------------------------------------------
    nexttile

    plot(timeSec, exampleTs, ...
        'LineWidth',1.1);

    hold on

    plot(timeSec, exampleRepl, ...
        'LineWidth',1.1);

    scatter( ...
        timeSec(exampleSpikeIdx), ...
        exampleRepl(exampleSpikeIdx), ...
        24, ...
        'filled');

    xlabel('Time (s)');
    ylabel('BOLD signal (a.u.)');

    legend( ...
        {'Original','Despiked','Edited samples'}, ...
        'Location','best');

    title(sprintf( ...
        ['Largest-edit example: voxel [%d %d %d], ' ...
         'max |change| = %.4g'], ...
        exampleVoxelSub(1), ...
        exampleVoxelSub(2), ...
        exampleVoxelSub(3), ...
        exampleMaxAbsChange));

    xlim([timeSec(1) timeSec(end)]);
    grid on


    % ------------------------------------------------------------
    % Figure title + output
    % ------------------------------------------------------------
    title(tl, sprintf( ...
        ['GM local edit QC: cut2 = %.3g | ' ...
         '%d samples edited (%.4f%% of GM voxel-time)'], ...
        cut2, ...
        nAlteredSamples, ...
        100 * nAlteredSamples / (nVox*T)));

    qcPng = fullfile(outputDir, ...
        [name, '_localGM_c2-', cutTag, '_QC.png']);

    exportgraphics(hFig, qcPng, ...
        'Resolution',180);

    close(hFig);

    qc.qc_plot_file = qcPng;

    fprintf('QC plot: %s\n', qcPng);

else

    qc.qc_plot_file = '';

end

%% Save output NIfTI, preserving original storage representation

scaledBack = (Y_out - offset) ./ scale;
outputData = cast(scaledBack, storedClass);

outNii = fullfile(outputDir, ...
    [name, '_despiked_localGM_c2-', cutTag]);

niftiwrite(outputData, outNii, funcInfo, "Compressed", true);

despiked_file = [outNii, '.nii.gz'];

fprintf('Output: %s\n', despiked_file);

end


%% Local helper: Preliminary spike squashing for fitting
%
% The AFNI source computes, at every time point, the median and MAD of a
% 9-point window. It then takes the median of those local MAD values,
% multiplies it by 6.789, and temporarily replaces values farther than that
% threshold from the local median. The modified series is used only to fit
% the smooth baseline.

function fitInput = cicada_new_prefit(ts, windowIdx)

windows = ts(windowIdx);

localMedian = median(windows, 2);
localMAD = median(abs(windows - localMedian), 2);

globalLocalMAD = median(localMAD);

fitInput = ts;

if ~isfinite(globalLocalMAD) || globalLocalMAD <= 0
    return;
end

preCut = 6.789 * globalLocalMAD;

preSpike = abs(ts - localMedian) > preCut;

fitInput(preSpike) = localMedian(preSpike);

end


%% Local helper: filename-safe numeric tag

function tag = cicada_numeric_tag(x)

tag = sprintf('%.6g', x);
tag = strrep(tag, '.', 'p');
tag = strrep(tag, '-', 'm');

end
