function [new_GM_probseg_file, new_WM_probseg_file, new_CSF_probseg_file, func_GM_file, func_WM_file, func_CSF_file, func_SUSC_file] = ...
    refineTPM(T1file, anatMaskFile, GMfile, WMfile, CSFfile, funcFile, outPrefix, maxIter)
% Refine TPMs (GM, WM, CSF) using T1 intensity + priors (Stage 1),
% then optionally refine TPMs in functional space with an added
% susceptibility class (Stage 2).
%
% Inputs:
%   T1file       = path to T1 NIfTI
%   anatMaskFile = path to anatomical mask NIfTI (binary or probabilistic)
%   GMfile       = path to GM TPM NIfTI
%   WMfile       = path to WM TPM NIfTI
%   CSFfile      = path to CSF TPM NIfTI
%   funcFile     = optional, path to 4D functional NIfTI (enables Stage 2)
%   outPrefix    = prefix for output files
%   maxIter      = optional, # iterations for Stage 1 (default 1)
%
% Outputs:
%   new_*_probseg_file = structurally refined TPMs (GM, WM, CSF) [Stage 1]
%   func_*_file        = functionally refined TPMs (GM, WM, CSF, SUSC) [Stage 2]
%
% Notes:
%  - Stage 1 is your original T1-only refinement.
%  - Stage 2 uses (i) the mean functional image, (ii) a graded susceptibility TPM
%    derived from inverse intensity, and (iii) the Stage-1 TPMs resampled to
%    functional space as priors. Then a 4-class GMM is fit on the functional mean.

if nargin < 7 || isempty(maxIter)
    maxIter = 1;
end

% === Stage 1: Structural refinement (unchanged logic) ====================
[GMfilepath, ~]  = fileparts(GMfile);
[WMfilepath, ~]  = fileparts(WMfile);
[CSFfilepath, ~] = fileparts(CSFfile);

% --- Load images with scaling ---
T1   = load_nii_scaled(T1file);
GM0  = load_nii_scaled(GMfile);
WM0  = load_nii_scaled(WMfile);
CSF0 = load_nii_scaled(CSFfile);
mask = load_nii_scaled(anatMaskFile) > 0; % logical mask

% Mask T1 and TPMs
T1   = T1  .* mask;
GM0  = GM0 .* mask;
WM0  = WM0 .* mask;
CSF0 = CSF0.* mask;

% Normalize priors
priorSum = GM0 + WM0 + CSF0 + eps;
GM0  = GM0  ./ priorSum;
WM0  = WM0  ./ priorSum;
CSF0 = CSF0 ./ priorSum;

% Initialize previous TPMs for convergence
prevGM = GM0; prevWM = WM0; prevCSF = CSF0;
tol = 1e-4;

% Vectorize masked data for GMM fitting
maskData = mask(:) > 0;
dataFit  = T1(maskData);

for it = 1:maxIter
    % --- Component stats on masked voxels ---
    [muGM, varGM]   = componentStats(dataFit, prevGM(maskData));
    [muWM, varWM]   = componentStats(dataFit, prevWM(maskData));
    [muCSF, varCSF] = componentStats(dataFit, prevCSF(maskData));

    startStruct = struct('mu', [muGM; muWM; muCSF], ...
                         'Sigma', reshape([varGM varWM varCSF], [1 1 3]), ...
                         'ComponentProportion', [1/3 1/3 1/3]);

    gmDist = fitgmdist(dataFit, 3, ...
                       'Start', startStruct, ...
                       'CovarianceType', 'diagonal', ...
                       'RegularizationValue', 1e-6, ...
                       'Options', statset('MaxIter',200));

    % Posteriors over the whole volume (unmasked)
    P = posterior(gmDist, T1(:));
    P = reshape(P, [size(T1,1), size(T1,2), size(T1,3), 3]);

    GM  = P(:,:,:,1) .* prevGM;
    WM  = P(:,:,:,2) .* prevWM;
    CSF = P(:,:,:,3) .* prevCSF;

    sumAll = GM + WM + CSF + eps;
    GM  = GM  ./ sumAll;
    WM  = WM  ./ sumAll;
    CSF = CSF ./ sumAll;

    % Convergence
    delta = mean(abs(GM(:)-prevGM(:))) + ...
            mean(abs(WM(:)-prevWM(:))) + ...
            mean(abs(CSF(:)-prevCSF(:)));

    if delta < tol
        fprintf('Stage 1 converged after %d iterations.\n', it);
        break;
    end

    prevGM = GM; prevWM = WM; prevCSF = CSF;
end

% Save Stage 1 outputs
infoT1 = niftiinfo(T1file);
infoT1.MultiplicativeScaling = 1;
infoT1.AdditiveOffset = 0;

new_GM_probseg_file  = fullfile(GMfilepath,  [outPrefix '_GM_probseg.nii.gz']);
new_WM_probseg_file  = fullfile(WMfilepath,  [outPrefix '_WM_probseg.nii.gz']);
new_CSF_probseg_file = fullfile(CSFfilepath, [outPrefix '_CSF_probseg.nii.gz']);

niftiwrite(single(GM),  strip_ext(fullfile(GMfilepath,  [outPrefix '_GM_probseg'])),  infoT1, 'Compressed', true);
niftiwrite(single(WM),  strip_ext(fullfile(WMfilepath,  [outPrefix '_WM_probseg'])),  infoT1, 'Compressed', true);
niftiwrite(single(CSF), strip_ext(fullfile(CSFfilepath, [outPrefix '_CSF_probseg'])), infoT1, 'Compressed', true);

% === Stage 2: Functional refinement (optional) ===========================
func_GM_file = ''; func_WM_file = ''; func_CSF_file = ''; func_SUSC_file = '';
if nargin < 8 || isempty(funcFile)
    fprintf('No functional file provided: skipping Stage 2 functional refinement.\n');
    return;
end

% Load functional (mean image) + mask to func space
funcInfo = niftiinfo(funcFile);
funcVol  = double(niftiread(funcInfo));
if ndims(funcVol) ~= 4
    error('funcFile should be 4D (x,y,z,t).');
end
funcMean = mean(funcVol, 4);

% resample structural-refined TPMs + anatomical mask to functional space if needed
GM_s  = maybe_resample_like(funcMean, funcInfo, GM,  infoT1);
WM_s  = maybe_resample_like(funcMean, funcInfo, WM,  infoT1);
CSF_s = maybe_resample_like(funcMean, funcInfo, CSF, infoT1);
mask_f = maybe_resample_like(funcMean, funcInfo, mask, niftiinfo(anatMaskFile)) > 0.5;

% Clean zeros; mask functional mean
funcMean = funcMean .* mask_f;

% --- Build graded susceptibility prior from inverse intensity -------------
% Exclude zeros/outside; robust percentile mapping (p5->1, p25->0, linear)
vals = funcMean(mask_f & funcMean > 0);
if isempty(vals)
    error('Functional mean has no positive voxels inside mask.');
end
p5  = prctile(vals, 5);
p25 = prctile(vals, 25);
den = max(p25 - p5, eps);
suscPrior = (p25 - funcMean) ./ den;         % linear map
suscPrior = min(max(suscPrior, 0), 1);       % clamp 0..1
suscPrior(~mask_f) = 0;

% optional light smoothing to avoid harsh edges
try
    suscPrior = imgaussfilt3(suscPrior, 0.5, 'FilterSize', 5);
catch
    % if Image Processing Toolbox absent, skip smoothing
end
suscPrior = min(max(suscPrior, 0), 1);

% --- Assemble 4-class priors and normalize --------------------------------
sum4 = GM_s + WM_s + CSF_s + suscPrior + eps;
GM_p   = GM_s   ./ sum4;
WM_p   = WM_s   ./ sum4;
CSF_p  = CSF_s  ./ sum4;
SUSC_p = suscPrior ./ sum4;

% --- Fit 4-class univariate GMM on funcMean (masked) ----------------------
maskDataF = mask_f(:) & (funcMean(:) > 0);
dataFitF  = funcMean(maskDataF);

% Component stats seeded from priors within func space
[muGMf, varGMf]     = componentStats(dataFitF, GM_p(maskDataF));
[muWMf, varWMf]     = componentStats(dataFitF, WM_p(maskDataF));
[muCSFf, varCSFf]   = componentStats(dataFitF, CSF_p(maskDataF));
[muSUSCf, varSUSCf] = componentStats(dataFitF, SUSC_p(maskDataF));

start4 = struct('mu', [muGMf; muWMf; muCSFf; muSUSCf], ...
                'Sigma', reshape([varGMf varWMf varCSFf varSUSCf], [1 1 4]), ...
                'ComponentProportion', [0.3 0.3 0.3 0.1]); % small prior mass to SUSC

gm4 = fitgmdist(dataFitF, 4, ...
                'Start', start4, ...
                'CovarianceType', 'diagonal', ...
                'RegularizationValue', 1e-6, ...
                'Options', statset('MaxIter',200));

% Posteriors over entire volume
P4 = nan(numel(funcMean), 4);
P4(maskDataF, :) = posterior(gm4, dataFitF);
P4(~maskDataF, :) = 0;
P4 = reshape(P4, [size(funcMean) 4]);

% Fuse with priors (Bayesian product), normalize to sum 1
GM_f   = P4(:,:,:,1) .* GM_p;
WM_f   = P4(:,:,:,2) .* WM_p;
CSF_f  = P4(:,:,:,3) .* CSF_p;
SUSC_f = P4(:,:,:,4) .* SUSC_p;

sumAllF = GM_f + WM_f + CSF_f + SUSC_f + eps;
GM_f   = GM_f   ./ sumAllF;
WM_f   = WM_f   ./ sumAllF;
CSF_f  = CSF_f  ./ sumAllF;
SUSC_f = SUSC_f ./ sumAllF;

% --- Save Stage 2 outputs -------------------------------------------------
funcOutDir = fileparts(funcFile);
funcInfo.MultiplicativeScaling = 1;
funcInfo.AdditiveOffset = 0;
funcInfo.ImageSize = funcInfo.ImageSize(1:3);
funcInfo.PixelDimensions = funcInfo.PixelDimensions(1:3);
funcInfo.Datatype = 'single';
funcInfo.BitsPerPixel = 16;

func_GM_file   = fullfile(funcOutDir, [outPrefix '_func_GM_probseg.nii.gz']);
func_WM_file   = fullfile(funcOutDir, [outPrefix '_func_WM_probseg.nii.gz']);
func_CSF_file  = fullfile(funcOutDir, [outPrefix '_func_CSF_probseg.nii.gz']);
func_SUSC_file = fullfile(funcOutDir, [outPrefix '_func_SUSC_probseg.nii.gz']);

niftiwrite(single(GM_f),   strip_ext(fullfile(funcOutDir, [outPrefix '_func_GM_probseg'])),   funcInfo, 'Compressed', true);
niftiwrite(single(WM_f),   strip_ext(fullfile(funcOutDir, [outPrefix '_func_WM_probseg'])),   funcInfo, 'Compressed', true);
niftiwrite(single(CSF_f),  strip_ext(fullfile(funcOutDir, [outPrefix '_func_CSF_probseg'])),  funcInfo, 'Compressed', true);
niftiwrite(single(SUSC_f), strip_ext(fullfile(funcOutDir, [outPrefix '_func_SUSC_probseg'])), funcInfo, 'Compressed', true);

end

%% --- Helper: compute mean/variance using top-weighted voxels -------------
function [muComp, varComp] = componentStats(data, prior)
    prior = prior(:);
    data  = data(:);
    if isempty(data) || isempty(prior)
        muComp = 0; varComp = 1;
        return;
    end
    prior = max(prior, 0);
    if all(prior==0)
        % Fallback to robust stats
        muComp = median(data);
        madv = mad(data,1) + 1e-4;
        varComp = max(madv.^2, 1e-4);
        return;
    end
    % take the top ~50% by prior weight (or at least 50 voxels)
    nTop = max(floor(0.5 * numel(prior)), 50);
    [~, idx] = sort(prior, 'descend');
    sel = idx(1:min(nTop, numel(idx)));
    w = prior(sel) + eps;
    d = data(sel);
    muComp  = sum(d .* w) / sum(w);
    varComp = sum(w .* (d - muComp).^2) / sum(w) + 1e-4;
end

%% --- Helper: load and scale NIfTI ---------------------------------------
function V = load_nii_scaled(filename)
    info = niftiinfo(filename);
    V = double(niftiread(info));
    if isfield(info, 'MultiplicativeScaling') && ~isempty(info.MultiplicativeScaling) && info.MultiplicativeScaling ~= 0
        V = V * info.MultiplicativeScaling;
    end
    if isfield(info, 'AdditiveOffset') && ~isempty(info.AdditiveOffset) && info.AdditiveOffset ~= 0
        V = V + info.AdditiveOffset;
    end
end

%% --- Helper: maybe resample source volume to match target geometry -------
function Vout = maybe_resample_like(targetVol, targetInfo, sourceVol, sourceInfo)
    % If size already matches, return as-is
    if isequal(size(targetVol,1), size(sourceVol,1)) && ...
       isequal(size(targetVol,2), size(sourceVol,2)) && ...
       isequal(size(targetVol,3), size(sourceVol,3))
        Vout = sourceVol;
        return;
    end
    % Simple trilinear resize (assumes roughly aligned axes)
    try
        Vout = imresize3(sourceVol, size(targetVol), 'linear');
    catch
        % Fallback if imresize3 not available
        [X,Y,Z] = ndgrid( linspace(1,size(sourceVol,1),size(targetVol,1)), ...
                          linspace(1,size(sourceVol,2),size(targetVol,2)), ...
                          linspace(1,size(sourceVol,3),size(targetVol,3)) );
        Vout = interpn(sourceVol, X, Y, Z, 'linear', 0);
    end
    % Clamp probabilities if it is TPM-like
    if max(Vout(:)) <= 1.5 && min(Vout(:)) >= -0.5
        Vout = min(max(Vout, 0), 1);
    end
end

%% --- Helper: strip extension (handles .nii and .nii.gz) ------------------
function pathNoExt = strip_ext(pathIn)
    % Removes .nii or .nii.gz so niftiwrite can append correctly
    if endsWith(pathIn, '.nii.gz')
        pathNoExt = extractBefore(pathIn, strlength(pathIn)-6);
    elseif endsWith(pathIn, '.nii')
        pathNoExt = extractBefore(pathIn, strlength(pathIn)-4);
    else
        pathNoExt = pathIn;
    end
end
