function [denoised_Edge_Edge_corr, denoised_FD_GM_corr, denoised_DVARS_GM_corr, ...
    denoised_Outbrain_Outbrain_corr, denoised_WMCSF_WMCSF_corr, denoised_CSF_CSF_corr, ...
    denoised_NotGM_NotGM_corr, denoised_GM_GM_corr, denoised_Suscept_Suscept_corr, ...
    denoised_GM_mean, denoised_mean_var_table, ...
    compare_Edge_Edge_corr, compare_FD_GM_corr, compare_DVARS_GM_corr, ...
    compare_Outbrain_Outbrain_corr, compare_WMCSF_WMCSF_corr, compare_CSF_CSF_corr, ...
    compare_NotGM_NotGM_corr, compare_GM_GM_corr, compare_Suscept_Suscept_corr, ...
    compare_GM_mean, compare_mean_var_table, ...
    orig_Edge_Edge_corr, orig_FD_GM_corr, orig_DVARS_GM_corr, ...
    orig_Outbrain_Outbrain_corr, orig_WMCSF_WMCSF_corr, orig_CSF_CSF_corr, ...
    orig_NotGM_NotGM_corr, orig_GM_GM_corr, orig_Suscept_Suscept_corr, ...
    orig_GM_mean, orig_mean_var_table] = CICADA_fileQC(cleaned_dir, denoised_file, compare_file, orig_file)
% An Inner Base script, called upon by script 3, to get relevant QC
% information
% Files should be in the cicada cleaned directory (for pathing purposes!)


if ~isfile(denoised_file)
    fprintf(['Cannot find denoised file at ', denoised_file, '\n'])
    return;
end

if ~isfile(compare_file)
    fprintf(['Cannot find compare file at ', compare_file, '\n'])
    return;
end

if ~isfile(orig_file)
    fprintf(['Cannot find orig file at ', orig_file, '\n'])
    return;
end

denoised_file_info = dir(denoised_file);
cd(cleaned_dir)

denoised_file_data = niftiread(denoised_file);
compare_file_data = niftiread(compare_file);
orig_file_data = niftiread(orig_file);

% confounds:
if exist([cleaned_dir, '/../confounds_timeseries.tsv'], 'file') ~= 0
    confound_place = [cleaned_dir, '/../confounds_timeseries.tsv'];
    % if it is a .tsv, then it should be tab delimited.
    allconfounds = readtable(confound_place, 'FileType', 'text', 'Delimiter', '\t');
elseif exist([cleaned_dir, '/../confounds_timeseries.csv'], 'file') ~= 0
    confound_place = [cleaned_dir, '/../confounds_timeseries.csv'];
    % matlab can natively read csv no problem. CSV is probably better.
    allconfounds = readtable(confound_place);
else
    fprintf('ERROR: Cannot find confounds_timeseries as either csv or tsv?')
    return;
end

confounds_dvars = table2array(allconfounds(:,{'dvars'}));
confounds_fd = table2array(allconfounds(:,{'framewise_displacement'}));

% Get GM prob mask and create GM mask and not GM mask for denoised
% and original files
GM_prob_file = [cleaned_dir, '/../region_masks/GM_prob.nii.gz']; % can reach into WM that directly neighbors GM, which often corresponds with smooth signal
WMCSF_prob_file = [cleaned_dir, '/../region_masks/Subepe_prob.nii.gz'];
CSF_prob_file = [cleaned_dir, '/../region_masks/CSF_prob.nii.gz'];
Edge_prob_file = [cleaned_dir, '/../region_masks/Edge_prob.nii.gz'];
Outbrain_prob_file = [cleaned_dir, '/../region_masks/OutbrainOnly_prob.nii.gz'];
Suscept_prob_file = [cleaned_dir, '/../region_masks/Susceptibility_prob.nii.gz'];
funcmask = [cleaned_dir, '/../funcmask.nii.gz'];
NotGM_prob_file = [cleaned_dir, '/../region_masks/NotGMorWM_prob.nii.gz']; % this is a better estimate if we calculated it, won't penalize (or give benefit to) strictly WM

fprintf(['Running QC calculation for ', denoised_file_info.name, '\n'])    
% Now, get data for NotGM, Edge, Outbrain, WMCSF, and CSF for the denoised
% data, compare data, and orig data
fprintf('Calculating Denoised Data\n')
% [denoised_GM, denoised_GM_mean, denoised_GM_mean_var] = getData_sameregion(denoised_file_data, funcmask, GM_prob_file);
% [denoised_GM_NotGM, denoised_NotGM_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, NotGM_prob_file);
% [denoised_GM_Edge, denoised_Edge_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Edge_prob_file);
% [denoised_GM_Outbrain, denoised_Outbrain_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Outbrain_prob_file);
% [denoised_GM_WMCSF, denoised_WMCSF_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, WMCSF_prob_file);
% [denoised_GM_CSF, denoised_CSF_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, CSF_prob_file);
% [denoised_GM_Suscept, denoised_Suscept_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Suscept_prob_file);

% Denoised Data in a given region
[denoised_GM, denoised_GM_mean, denoised_GM_mean_var] = getData_sameregion(denoised_file_data, funcmask, GM_prob_file);
[denoised_NotGM, ~, denoised_NotGM_mean_var] = getData_sameregion(denoised_file_data, funcmask, NotGM_prob_file);
[denoised_Edge, ~, denoised_Edge_mean_var] = getData_sameregion(denoised_file_data, funcmask, Edge_prob_file);
[denoised_Outbrain, ~, denoised_Outbrain_mean_var] = getData_sameregion(denoised_file_data, funcmask, Outbrain_prob_file);
[denoised_WMCSF, ~, denoised_WMCSF_mean_var] = getData_sameregion(denoised_file_data, funcmask, WMCSF_prob_file);
[denoised_CSF, ~, denoised_CSF_mean_var] = getData_sameregion(denoised_file_data, funcmask, CSF_prob_file);
[denoised_Suscept, ~, denoised_Suscept_mean_var] = getData_sameregion(denoised_file_data, funcmask, Suscept_prob_file);

% Compare Data
[compare_GM, compare_GM_mean, compare_GM_mean_var] = getData_sameregion(compare_file_data, funcmask, GM_prob_file);
[compare_NotGM, ~, compare_NotGM_mean_var] = getData_sameregion(compare_file_data, funcmask, NotGM_prob_file);
[compare_Edge, ~, compare_Edge_mean_var] = getData_sameregion(compare_file_data, funcmask, Edge_prob_file);
[compare_Outbrain, ~, compare_Outbrain_mean_var] = getData_sameregion(compare_file_data, funcmask, Outbrain_prob_file);
[compare_WMCSF, ~, compare_WMCSF_mean_var] = getData_sameregion(compare_file_data, funcmask, WMCSF_prob_file);
[compare_CSF, ~, compare_CSF_mean_var] = getData_sameregion(compare_file_data, funcmask, CSF_prob_file);
[compare_Suscept, ~, compare_Suscept_mean_var] = getData_sameregion(compare_file_data, funcmask, Suscept_prob_file);

% Original Data
[orig_GM, orig_GM_mean, orig_GM_mean_var] = getData_sameregion(orig_file_data, funcmask, GM_prob_file);
[orig_NotGM, ~, orig_NotGM_mean_var] = getData_sameregion(orig_file_data, funcmask, NotGM_prob_file);
[orig_Edge, ~, orig_Edge_mean_var] = getData_sameregion(orig_file_data, funcmask, Edge_prob_file);
[orig_Outbrain, ~, orig_Outbrain_mean_var] = getData_sameregion(orig_file_data, funcmask, Outbrain_prob_file);
[orig_WMCSF, ~, orig_WMCSF_mean_var] = getData_sameregion(orig_file_data, funcmask, WMCSF_prob_file);
[orig_CSF, ~, orig_CSF_mean_var] = getData_sameregion(orig_file_data, funcmask, CSF_prob_file);
[orig_Suscept, ~, orig_Suscept_mean_var] = getData_sameregion(orig_file_data, funcmask, Suscept_prob_file);


% Now compute relevant correlations, and some others that we don't use, but
% one could make use of if they desired
fprintf('Calculating Relevant Correlations\n')

% create permutations for regions, subtract a little bit because sometimes
% something funky happens and the sizes do not fully align with Edge
% Voxels. Very close (probably related to thresholding values and minor 
% smoothing meaning some voxels that are all 0s for one, are not for the other.)

% For GM compared to other region
GM_randperm = randperm(size(denoised_GM, 1));
% GM_Edge_randperm = randperm(size(denoised_GM_Edge, 1));
% GM_Outbrain_randperm = randperm(size(denoised_GM_Outbrain, 1));
% GM_WMCSF_randperm = randperm(size(denoised_GM_WMCSF, 1));
% GM_CSF_randperm = randperm(size(denoised_GM_CSF, 1));
% GM_Suscept_randperm = randperm(size(denoised_GM_Suscept, 1));
% GM_NotGM_randperm = randperm(size(denoised_GM_NotGM, 1));

% For other region compared to GM
% NotGM_GM_randperm = randperm(size(denoised_NotGM_GM, 1));
% Outbrain_GM_randperm = randperm(size(denoised_Outbrain_GM, 1));
% Edge_GM_randperm = randperm(size(denoised_Edge_GM, 1));
% WMCSF_GM_randperm = randperm(size(denoised_WMCSF_GM, 1));
% CSF_GM_randperm = randperm(size(denoised_CSF_GM, 1));
% Suscept_GM_randperm = randperm(size(denoised_Suscept_GM, 1));

% for region compared to itself:
NotGM_randperm = randperm(size(denoised_NotGM, 1));
Outbrain_randperm = randperm(size(denoised_Outbrain, 1));
Edge_randperm = randperm(size(denoised_Edge, 1));
WMCSF_randperm = randperm(size(denoised_WMCSF, 1));
CSF_randperm = randperm(size(denoised_CSF, 1));
Suscept_randperm = randperm(size(denoised_Suscept, 1));

% Now get all the relevant correlations.
%
% IMPORTANT: CICADA, comparison, and original data must use the same
% sampled voxels / voxel pairs so that subsequent paired statistical tests
% are truly paired. Undefined correlations are removed jointly across all
% three datasets.
%
% perms should generally be at least 1000, but can go lower if a region
% contains fewer usable voxels/pairs.
perms = 10000;


% Edge
[denoised_Edge_Edge_corr, ...
    compare_Edge_Edge_corr, ...
    orig_Edge_Edge_corr] = ...
    createHistDataPaired( ...
    denoised_Edge, compare_Edge, orig_Edge, ...
    Edge_randperm, perms);


% FD vs GM
[denoised_FD_GM_corr, ...
    compare_FD_GM_corr, ...
    orig_FD_GM_corr] = ...
    createHistConf1DDataPaired( ...
    denoised_GM, compare_GM, orig_GM, ...
    confounds_fd, GM_randperm, perms);


% DVARS vs GM
[denoised_DVARS_GM_corr, ...
    compare_DVARS_GM_corr, ...
    orig_DVARS_GM_corr] = ...
    createHistConf1DDataPaired( ...
    denoised_GM, compare_GM, orig_GM, ...
    confounds_dvars, GM_randperm, perms);


% Outbrain
[denoised_Outbrain_Outbrain_corr, ...
    compare_Outbrain_Outbrain_corr, ...
    orig_Outbrain_Outbrain_corr] = ...
    createHistDataPaired( ...
    denoised_Outbrain, compare_Outbrain, orig_Outbrain, ...
    Outbrain_randperm, perms);


% WM / CSF
[denoised_WMCSF_WMCSF_corr, ...
    compare_WMCSF_WMCSF_corr, ...
    orig_WMCSF_WMCSF_corr] = ...
    createHistDataPaired( ...
    denoised_WMCSF, compare_WMCSF, orig_WMCSF, ...
    WMCSF_randperm, perms);


% CSF
[denoised_CSF_CSF_corr, ...
    compare_CSF_CSF_corr, ...
    orig_CSF_CSF_corr] = ...
    createHistDataPaired( ...
    denoised_CSF, compare_CSF, orig_CSF, ...
    CSF_randperm, perms);


% Susceptibility-prone region
[denoised_Suscept_Suscept_corr, ...
    compare_Suscept_Suscept_corr, ...
    orig_Suscept_Suscept_corr] = ...
    createHistDataPaired( ...
    denoised_Suscept, compare_Suscept, orig_Suscept, ...
    Suscept_randperm, perms);


% Non-GM
[denoised_NotGM_NotGM_corr, ...
    compare_NotGM_NotGM_corr, ...
    orig_NotGM_NotGM_corr] = ...
    createHistDataPaired( ...
    denoised_NotGM, compare_NotGM, orig_NotGM, ...
    NotGM_randperm, perms);


% GM
[denoised_GM_GM_corr, ...
    compare_GM_GM_corr, ...
    orig_GM_GM_corr] = ...
    createHistDataPaired( ...
    denoised_GM, compare_GM, orig_GM, ...
    GM_randperm, perms);

%Edge_GM_corr = createHistData(denoised_GM_Edge, denoised_Edge_GM, GM_Edge_randperm, Edge_GM_randperm, perms);
%Outbrain_GM_corr = createHistData(denoised_GM_Outbrain, denoised_Outbrain_GM, GM_Outbrain_randperm, Outbrain_GM_randperm, perms);
%WMCSF_GM_corr = createHistData(denoised_GM_WMCSF, denoised_WMCSF_GM, GM_WMCSF_randperm, WMCSF_GM_randperm, perms);
%CSF_GM_corr = createHistData(denoised_GM_CSF, denoised_CSF_GM, GM_CSF_randperm, CSF_GM_randperm, perms);
%NotGM_GM_corr = createHistData(denoised_GM_NotGM, denoised_NotGM_GM, GM_NotGM_randperm, NotGM_GM_randperm, perms);
%Suscept_GM_corr = createHistData(denoised_GM_Suscept, denoised_Suscept_GM, GM_Suscept_randperm, Suscept_GM_randperm, perms);

% denoised
denoised_mean_var_array = [denoised_GM_mean_var, denoised_NotGM_mean_var, denoised_Edge_mean_var, ...
    denoised_Outbrain_mean_var, denoised_WMCSF_mean_var, denoised_CSF_mean_var, denoised_Suscept_mean_var];
mean_var_labels = ["GM_mean_var", "NotGM_mean_var", "Edge_mean_var", "Outbrain_mean_var", ...
    "WMCSF_mean_var", "CSF_mean_var", "Suscept_mean_var"];
denoised_mean_var_table = array2table(denoised_mean_var_array, 'VariableNames', mean_var_labels);

% compare
compare_mean_var_array = [compare_GM_mean_var, compare_NotGM_mean_var, compare_Edge_mean_var, ...
    compare_Outbrain_mean_var, compare_WMCSF_mean_var, compare_CSF_mean_var, compare_Suscept_mean_var];
mean_var_labels = ["GM_mean_var", "NotGM_mean_var", "Edge_mean_var", "Outbrain_mean_var", ...
    "WMCSF_mean_var", "CSF_mean_var", "Suscept_mean_var"];
compare_mean_var_table = array2table(compare_mean_var_array, 'VariableNames', mean_var_labels);

% orig
orig_mean_var_array = [orig_GM_mean_var, orig_NotGM_mean_var, orig_Edge_mean_var, ...
    orig_Outbrain_mean_var, orig_WMCSF_mean_var, orig_CSF_mean_var, orig_Suscept_mean_var];
mean_var_labels = ["GM_mean_var", "NotGM_mean_var", "Edge_mean_var", "Outbrain_mean_var", ...
    "WMCSF_mean_var", "CSF_mean_var", "Suscept_mean_var"];
orig_mean_var_table = array2table(orig_mean_var_array, 'VariableNames', mean_var_labels);

end

%%

function [region, region_signal_mean, region_mean_var] = getData_sameregion(funcfile_data, funcmask, region_prob_file)
% read in the niftis
region_prob_data = niftiread(region_prob_file);
func_data = funcfile_data;
funcmask_data = niftiread(funcmask);

region_mask = logical(region_prob_data > 0.67 & funcmask_data == 1);

% Select the timeseries for GM and region
region = zeros(sum(region_mask(:)), size(func_data, 4));
for j = 1:size(func_data, 4)
    curr_func = func_data(:,:,:,j);
    region(:,j) = curr_func(region_mask);
end

region_mean_var = mean(var(region, 0, 2), "omitnan");
region_signal_mean = mean(region(region(:,1) ~= 0, :));
region_signal_mean = region_signal_mean - mean(region_signal_mean);
end

function [denoised_corr, compare_corr, orig_corr] = ...
    createHistDataPaired(denoised_region, compare_region, orig_region, ...
    region_randperm, perms)
% Sample the same voxel pairs from denoised, comparison, and original data.
%
% This preserves observation-level pairing for downstream paired tests.
% Rather than forming a very large all-by-all correlation matrix and then
% independently sampling it for each dataset, randomly sample the desired
% voxel pairs first and calculate only those correlations.

n_denoised = size(denoised_region, 1);
n_compare = size(compare_region, 1);
n_orig = size(orig_region, 1);

% These datasets are extracted using the same spatial region mask, so
% paired QC requires identical voxel counts/order.
if n_denoised ~= n_compare || n_denoised ~= n_orig

    error('CICADA:QCVoxelCountMismatch', ...
        ['Paired QC requires the denoised, comparison, and original ' ...
         'datasets to contain the same region voxels. Counts were ' ...
         '%d, %d, and %d.'], ...
        n_denoised, n_compare, n_orig);

end

nvox = min([n_denoised, numel(region_randperm), perms]);

if nvox < 2

    error('CICADA:QCInsufficientVoxels', ...
        'At least two region voxels are required for voxel-pair QC.');

end

% Use the same randomized voxel subset for all three datasets.
voxel_idx = region_randperm(1:nvox);

% Number of unique unordered voxel pairs.
num_pairs = nvox * (nvox - 1) / 2;

% Preserve the previous intent of sampling at most "perms" correlations.
n_sample = min(perms, num_pairs);

% Uniformly sample unique unordered pairs without first creating the full
% nvox-by-nvox correlation matrix.
pair_rank = randperm(num_pairs, n_sample);

% Convert pair ranks 1..nchoosek(nvox,2) into row-based unordered pairs:
% (2,1), (3,1), (3,2), (4,1), ...
pair_i = ceil((1 + sqrt(1 + 8 .* pair_rank)) ./ 2);
pair_j = pair_rank - ((pair_i - 1) .* (pair_i - 2) ./ 2);

voxel_i = voxel_idx(pair_i);
voxel_j = voxel_idx(pair_j);

% Calculate exactly the same voxel-pair correlations in all three datasets.
denoised_corr = local_rowwise_corr( ...
    denoised_region(voxel_i, :), ...
    denoised_region(voxel_j, :));

compare_corr = local_rowwise_corr( ...
    compare_region(voxel_i, :), ...
    compare_region(voxel_j, :));

orig_corr = local_rowwise_corr( ...
    orig_region(voxel_i, :), ...
    orig_region(voxel_j, :));

% A pair must be defined in all three datasets to remain in paired QC.
valid = ...
    isfinite(denoised_corr) & ...
    isfinite(compare_corr) & ...
    isfinite(orig_corr);

denoised_corr = denoised_corr(valid);
compare_corr = compare_corr(valid);
orig_corr = orig_corr(valid);

if isempty(denoised_corr)

    error('CICADA:QCNoValidCorrelations', ...
        ['No jointly valid voxel-pair correlations remained for this ' ...
         'QC region.']);

end

end


function [denoised_corr, compare_corr, orig_corr] = ...
    createHistConf1DDataPaired(denoised_GM, compare_GM, orig_GM, ...
    conf, GM_randperm, perms)
% Correlate FD/DVARS-like 1-D confounds with the same GM voxels in all
% datasets, preserving voxel-level pairing.
%
% Confounds are calculated in the first derivative / absolute-value domain
% (e.g., FD and DVARS), matching the historical CICADA QC calculation.

n_denoised = size(denoised_GM, 1);
n_compare = size(compare_GM, 1);
n_orig = size(orig_GM, 1);

if n_denoised ~= n_compare || n_denoised ~= n_orig

    error('CICADA:QCVoxelCountMismatch', ...
        ['Paired QC requires the denoised, comparison, and original ' ...
         'datasets to contain the same GM voxels. Counts were ' ...
         '%d, %d, and %d.'], ...
        n_denoised, n_compare, n_orig);

end

nvox = min([n_denoised, numel(GM_randperm), perms]);

if nvox < 1

    error('CICADA:QCInsufficientVoxels', ...
        'At least one GM voxel is required for confound QC.');

end

voxel_idx = GM_randperm(1:nvox);

% Same confound and same voxels for every denoising method.
denoised_corr = corr( ...
    conf(2:end), ...
    abs(diff(denoised_GM(voxel_idx, :)', 1, 1)))';

compare_corr = corr( ...
    conf(2:end), ...
    abs(diff(compare_GM(voxel_idx, :)', 1, 1)))';

orig_corr = corr( ...
    conf(2:end), ...
    abs(diff(orig_GM(voxel_idx, :)', 1, 1)))';

% Remove undefined correlations jointly so that observation pairing is
% retained for signrank and all downstream summaries.
valid = ...
    isfinite(denoised_corr) & ...
    isfinite(compare_corr) & ...
    isfinite(orig_corr);

denoised_corr = denoised_corr(valid);
compare_corr = compare_corr(valid);
orig_corr = orig_corr(valid);

if isempty(denoised_corr)

    error('CICADA:QCNoValidCorrelations', ...
        ['No jointly valid GM-confound correlations remained for this ' ...
         'QC metric.']);

end

end


function r = local_rowwise_corr(data_a, data_b)
% Pearson correlation between corresponding rows of two matrices.
%
% This avoids constructing a full voxel-by-voxel correlation matrix.

data_a = double(data_a);
data_b = double(data_b);

data_a = data_a - mean(data_a, 2);
data_b = data_b - mean(data_b, 2);

numerator = sum(data_a .* data_b, 2);

denominator = sqrt( ...
    sum(data_a .^ 2, 2) .* ...
    sum(data_b .^ 2, 2));

r = numerator ./ denominator;

% Constant/invalid time series have undefined correlation.
r(denominator == 0) = NaN;

end


