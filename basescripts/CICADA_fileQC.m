function [Edge_GM_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, NotGM_GM_corr, GM_GM_autocorr, GM_mean] = CICADA_fileQC(denoised_file, orig_file)
% An Inner Base script, called upon by script 3, to get relevant QC
% information
% denoised file and MUST be in the cicada cleaned directory (for pathing purposes!)
% we need orig file in order to do GM GM autocorrelation


if ~isfile(denoised_file)
    fprintf(['Cannot find denoised file at ', denoised_file, '\n'])
    return;
end

if ~isfile(orig_file)
    fprintf(['Cannot find orig file at ', orig_file, '\n'])
    return;
end

denoised_file_info = dir(denoised_file);
cleaned_dir = denoised_file_info.folder;
cd(cleaned_dir)

denoised_file_data = niftiread(denoised_file);
orig_file_data = niftiread(orig_file);

% confounds: So, overall you want/need the 6 motion parameters, dvars,
% framewise displacement, csf, white_matter, and global signal
confound_place = [cleaned_dir, '/../confounds_timeseries.csv'];
allconfounds = readtable(confound_place);
confounds_dvars = table2array(allconfounds(:,{'dvars'}));
confounds_fd = table2array(allconfounds(:,{'framewise_displacement'}));

% Get GM prob mask and create GM mask and not GM mask for denoised
% and original files
GM_prob_file = [cleaned_dir, '/../region_masks/GM_prob.nii.gz'];
WMCSF_prob_file = [cleaned_dir, '/../region_masks/WMandCSF_mask.nii.gz'];
CSF_prob_file = [cleaned_dir, '/../region_masks/CSF_prob.nii.gz'];
Edge_prob_file = [cleaned_dir, '/../region_masks/Edge_prob.nii.gz'];
Outbrain_prob_file = [cleaned_dir, '/../region_masks/OutbrainOnly_prob.nii.gz'];
funcmask = [cleaned_dir, '/../funcmask.nii.gz'];
if isfile([cleaned_dir, '/../region_masks/Outbrain_prob.nii.gz']) == 1
    NotGM_prob_file = [cleaned_dir, '/../region_masks/Outbrain_prob.nii.gz']; % this is a better estimate if we calculated it
else
    NotGM_prob_file = [cleaned_dir, '/../region_masks/NotGMorWM_prob.nii.gz']; % Makes sure not including anything that could be signal, removes WM too
end

fprintf(['Running QC calculation for ', denoised_file_info.name, '\n'])    

% Now, get data for NotGM, Edge, Outbrain, WMCSF, and CSF for the denoised
% data and orig data
fprintf('Calculating Denoised Data\n')
[denoised_GM, denoised_NotGM, ~, denoised_GM_mean, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, NotGM_prob_file);
[~, denoised_Edge, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Edge_prob_file);
[~, denoised_Outbrain, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Outbrain_prob_file);
[~, denoised_WMCSF, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, WMCSF_prob_file);
[~, denoised_CSF, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, CSF_prob_file);

% and then grab just GM for orig data to do auto corr 
[orig_GM, ~, ~, ~, ~, ~] = getData(orig_file_data, funcmask, GM_prob_file, NotGM_prob_file);

% Now compute relevant correlations, and some others that we don't use, but
% one could make use of if they desired
fprintf('Calculating Relevant Correlations\n')

% create permutations for regions, subtract a little bit because sometimes
% something funky happens and the sizes do not fully align with Edge
% Voxels. Very close (probably related to thresholding values and minor 
% smoothing meaning some voxels that are all 0s for one, are not for the other.)
GM_randperm = randperm(size(denoised_GM, 1));
NotGM_randperm = randperm(size(denoised_NotGM, 1));
Outbrain_randperm = randperm(size(denoised_Outbrain, 1));
Edge_randperm = randperm(size(denoised_Edge, 1));
WMCSF_randperm = randperm(size(denoised_WMCSF, 1));
CSF_randperm = randperm(size(denoised_CSF, 1));

% Now get all the relevant correlations:
Edge_GM_corr = createHistData(denoised_GM, denoised_Edge, GM_randperm, Edge_randperm, 500);
FD_GM_corr = createHistConf1DData(denoised_GM, confounds_fd, GM_randperm, 10000);
DVARS_GM_corr = createHistConf1DData(denoised_GM, confounds_dvars, GM_randperm, 10000);
Outbrain_GM_corr = createHistData(denoised_GM, denoised_Outbrain, GM_randperm, Outbrain_randperm, 500);
WMCSF_GM_corr = createHistData(denoised_GM, denoised_WMCSF, GM_randperm, WMCSF_randperm, 500);
CSF_GM_corr = createHistData(denoised_GM, denoised_CSF, GM_randperm, CSF_randperm, 500);
NotGM_GM_corr = createHistData(denoised_GM, denoised_NotGM, GM_randperm, NotGM_randperm, 500);
GM_GM_autocorr = createHistData(denoised_GM, orig_GM, GM_randperm, GM_randperm, 500); % if/when you feed it orig_file, then this will all just be 1

GM_mean = denoised_GM_mean;

end

%%
function [GM, region, Global_Signal, GM_Signal_mean, GM_var, region_var] = getData(funcfile_data, funcmask, GM_prob_file, region_prob_file)
% read in the niftis
region_prob_data = niftiread(region_prob_file);
gm_prob_data = niftiread(GM_prob_file);
func_data = funcfile_data;
funcmask_data = niftiread(funcmask);

gm_mask = logical(gm_prob_data > 0.75 & funcmask_data == 1);
region_mask = logical(region_prob_data > 0.75 & funcmask_data == 1);
func_mask = logical(funcmask_data == 1);

% Select the timeseries for GM and region and global signal
Global_Signal = zeros(sum(func_mask(:)), size(func_data, 4));
GM = zeros(sum(gm_mask(:)), size(func_data, 4));
region = zeros(sum(region_mask(:)), size(func_data, 4));
for j = 1:size(func_data, 4)
    curr_func = func_data(:,:,:,j);
    GM(:,j) = curr_func(gm_mask);
    region(:,j) = curr_func(region_mask);
    Global_Signal(:,j) = curr_func(func_mask);
end

%GM = GM(GM(:,1) ~= 0, :);
%region = region(region(:, 1) ~= 0, :);
%Global_Signal = global_signal(global_signal(:, 1) ~= 0, :);

% need to make sure size of GM and region are the same
if size(GM) ~= size(region)
    fprintf('Sizes of GM vs Region are different. This should not happen if made from the same mask?\n')
    return;
end

GM_var = var(GM, 0, 2); % variance of difference in time for GM)
region_var = var(region, 0, 2);

%Global_Signal_mean = mean(global_signal(global_signal(:, 1) ~= 0, :));
%Global_Signal_mean = Global_Signal_mean - mean(Global_Signal_mean);
GM_Signal_mean = mean(GM(GM(:,1) ~= 0, :));
GM_Signal_mean = GM_Signal_mean - mean(GM_Signal_mean);
end

function compare_orig_corr = createHistData(compare_GM, orig_region, GM_randperm, region_randperm, perms)
% Should center close to 0 (not globally correlated) if it's
% denoised (e.g., denoised lowered noise)
compare_orig_corr = corr(compare_GM(GM_randperm(1:perms), :)', orig_region(region_randperm(1:perms), :)');
compare_orig_corr = compare_orig_corr(tril(compare_orig_corr, -1) ~= 0);

% Should center above 0 (e.g., denoise did not remove true signal)
% GM_corr = corr(compare_GM(GM_randperm(1:perms), :)', compare_GM(GM_randperm(1:perms), :)');
% GM_GM_corr = GM_corr(tril(GM_corr, -1) ~= 0);
end

function [GM_conf1D_corr] = createHistConf1DData(GM, conf, GM_randperm, perms)
% This is for confounds that are calculated in the 1st derivative
% and absolute value (e.g., dvars, fd)
% corr conf vs GM diff abs: Should be reduced if denoised
% consider detrend() for GM, since that is what we do for script 2 to focus
% more on spikes, less on drift.
GM_conf1D_corr = corr(conf(2:end), abs(diff(GM(GM_randperm(1:perms), :)')))';
end