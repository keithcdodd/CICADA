function [dataDespiked, madMap, spikeCounts, logSummary] = despike_fMRI_derivative(data4D, varargin)
% despike_fMRI_derivative - Robust derivative-based despiking of 4D fMRI data
%
% Syntax:
%   [dataDespiked, madMap, spikeCounts, logSummary] = despike_fMRI_derivative(data4D, 'Param', Value, ...)
%
% Inputs:
%   data4D     - 4D fMRI data array (X x Y x Z x T)
%
% Optional Parameters (Name-Value pairs):
%   'ZThreshold'      - z-score threshold for spike detection (default: 6)
%   'Scale'           - softness scale factor for tanh attenuation (default: 2)
%   'Mask'            - 3D logical mask to restrict voxels despiked (default: all voxels)
%   'Verbose'         - logical, print progress info (default: true)
%   'SaveLogPath'     - full path to save despiking summary text file (default: '')
%
% Outputs:
%   dataDespiked - despiked 4D fMRI data array, same size as input
%   madMap       - 3D map of MAD-derived robust std dev of derivative for each voxel
%   spikeCounts  - 3D map of number of spikes detected per voxel
%   logSummary   - struct with summary statistics:
%                  .totalSpikes, .pctVoxelsWithSpikes, .meanMAD, .medianMAD, .maxMAD
%
% Description:
%   Computes temporal derivative per voxel, calculates robust z-scores
%   (median, scaled MAD), and softly attenuates spikes above threshold
%   with tanh-based compression. Minimal distortion otherwise.
%
% Example usage:
%   [dataDespiked, madMap, spikeCounts, logSummary] = despike_fMRI_derivative(data, ...
%       'ZThreshold', 6, 'Mask', gmMask, 'SaveLogPath', 'sub-01/despike_log.txt');
%
% Author: ChatGPT, adapted for your pipeline
% Date: 2025-07-06

% Parse inputs
p = inputParser;
addParameter(p, 'ZThreshold', 6, @(x) isnumeric(x) && isscalar(x) && (x>0));
addParameter(p, 'Scale', 2, @(x) isnumeric(x) && isscalar(x) && (x>0));
addParameter(p, 'Mask', [], @(x) (islogical(x) && ndims(x) == 3) || isempty(x));
addParameter(p, 'Verbose', true, @(x) islogical(x));
addParameter(p, 'SaveLogPath', '', @(x) ischar(x) || isstring(x));
parse(p, varargin{:});

zThresh = p.Results.ZThreshold;
scale = p.Results.Scale;
mask = p.Results.Mask;
verbose = p.Results.Verbose;
saveLogPath = p.Results.SaveLogPath;

if verbose
    fprintf('Starting despiking with z-threshold = %.2f, scale = %.2f\n', zThresh, scale);
end

[X, Y, Z, T] = size(data4D);
if T < 3
    error('Data time dimension too short to despike.');
end

% Default mask: all voxels
if isempty(mask)
    mask = true(X, Y, Z);
elseif ~isequal(size(mask), [X Y Z])
    error('Mask size must match spatial dimensions of data.');
end

dataDespiked = zeros(size(data4D), 'like', data4D);
madMap = zeros(X, Y, Z);
spikeCounts = zeros(X, Y, Z);

% Flatten spatial dims for speed (voxels x time)
data2D = reshape(data4D, [], T);
maskVec = reshape(mask, [], 1);
Nvox = sum(maskVec);

dataDespiked2D = zeros(size(data2D), 'like', data2D);

voxelCounter = 0;
totalSpikes = 0;
madVals = zeros(Nvox,1);

for v = 1:size(data2D,1)
    if ~maskVec(v)
        % Outside mask: copy original data
        dataDespiked2D(v, :) = data2D(v, :);
        continue;
    end
    
    voxelCounter = voxelCounter + 1;
    ts = double(data2D(v, :));  % numeric stability
    
    if std(ts) < eps
        dataDespiked2D(v, :) = ts;
        continue;
    end
    
    dts = diff(ts);
    
    meddts = median(dts);
    maddts_raw = mad(dts, 1);
    maddts = maddts_raw * 1.4826;  % normalized MAD ~ std dev
    
    if maddts < eps
        zscores = zeros(size(dts));
    else
        zscores = (dts - meddts) / maddts;
    end
    
    spikes = abs(zscores) > zThresh;
    nSpikes = sum(spikes);
    spikeCounts(v) = nSpikes;
    totalSpikes = totalSpikes + nSpikes;
    madVals(voxelCounter) = maddts;
    madMap(v) = maddts;
    
    % Soft threshold spikes with tanh scaling
    z_adj = zscores;
    z_adj(spikes) = sign(zscores(spikes)) .* zThresh .* tanh(abs(zscores(spikes))/scale);
    
    dts_adj = z_adj * maddts + meddts;
    
    ts_despiked = [ts(1), ts(1) + cumsum(dts_adj)];
    
    dataDespiked2D(v, :) = ts_despiked;
    
    if verbose && mod(voxelCounter, 10000) == 0
        fprintf('Processed %d/%d voxels (masked)\n', voxelCounter, Nvox);
    end
end

dataDespiked = reshape(dataDespiked2D, X, Y, Z, T);

pctVoxelsWithSpikes = 100 * nnz(spikeCounts) / Nvox;

logSummary.totalSpikes = totalSpikes;
logSummary.pctVoxelsWithSpikes = pctVoxelsWithSpikes;
logSummary.meanMAD = mean(madVals);
logSummary.medianMAD = median(madVals);
logSummary.maxMAD = max(madVals);

if verbose
    fprintf('Despiking complete.\n');
    fprintf('Total spikes detected: %d\n', totalSpikes);
    fprintf('Percent of voxels with spikes: %.2f%%\n', pctVoxelsWithSpikes);
    fprintf('MAD of derivative across voxels (mean/median/max): %.4f / %.4f / %.4f\n', ...
        logSummary.meanMAD, logSummary.medianMAD, logSummary.maxMAD);
end

% Save log summary text file if requested
if ~isempty(saveLogPath)
    try
        fid = fopen(saveLogPath, 'w');
        if fid == -1
            warning('Could not open file to save despiking log: %s', saveLogPath);
        else
            fprintf(fid, 'Despiking Summary Log\n');
            fprintf(fid, '---------------------\n');
            fprintf(fid, 'Date: %s\n', datestr(now));
            fprintf(fid, 'Z-Threshold: %.2f\n', zThresh);
            fprintf(fid, 'Softness Scale: %.2f\n', scale);
            fprintf(fid, 'Total voxels despiked: %d\n', Nvox);
            fprintf(fid, 'Total spikes detected: %d\n', totalSpikes);
            fprintf(fid, 'Percent voxels with spikes: %.2f%%\n', pctVoxelsWithSpikes);
            fprintf(fid, 'MAD derivative (mean/median/max): %.4f / %.4f / %.4f\n', ...
                logSummary.meanMAD, logSummary.medianMAD, logSummary.maxMAD);
            fclose(fid);
            if verbose
                fprintf('Despiking log saved to %s\n', saveLogPath);
            end
        end
    catch ME
        warning('Error saving despiking log: %s', ME.message);
    end
end

% Call like: saveFile = fullfile(subjectFolder, 'despike_summary.txt');
% but I like full brain bask instead of gmMask

%[dataDespiked, madMap, spikeCounts, logSummary] = despike_fMRI_derivative(data4D, ...
%    'ZThreshold', 6, 'Scale', 2, 'Mask', gmMask, 'SaveLogPath', saveFile, 'Verbose', true);


end
