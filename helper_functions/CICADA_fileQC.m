function [Edge_Edge_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_Outbrain_corr, WMCSF_WMCSF_corr, CSF_CSF_corr, NotGM_NotGM_corr, GM_GM_corr, Suscept_Suscept_corr, GM_mean] = CICADA_fileQC(denoised_file, orig_file)
% function [Edge_GM_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, NotGM_GM_corr, GM_GM_corr, GM_mean] = CICADA_fileQC(denoised_file, orig_file)
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

% consider meanfiltering in 3D the data, just to fill in potential holes
%smoothed_denoised_file_data = zeros(size(denoised_file_data));
%for idx = 1:size(denoised_file_data, 4)
%    smoothed_denoised_file_data(:,:,:,idx) = imgaussfilt3(denoised_file_data(:,:,:,idx),1); % small gauss filter
%end
%denoised_file_data = smoothed_denoised_file_data; % helps keep regions smoother, but not bleeding too much into other regions

% confounds: So, overall you want/need the 6 motion parameters, dvars,
% framewise displacement, csf, white_matter, and global signal
confound_place = [cleaned_dir, '/../confounds_timeseries.csv'];
allconfounds = readtable(confound_place);
confounds_dvars = table2array(allconfounds(:,{'dvars'}));
confounds_fd = table2array(allconfounds(:,{'framewise_displacement'}));

% Get GM prob mask and create GM mask and not GM mask for denoised
% and original files
GM_prob_file = [cleaned_dir, '/../region_masks/GM_extended_prob.nii.gz']; % reaches further into other regions
WMCSF_prob_file = [cleaned_dir, '/../region_masks/Subepe_prob.nii.gz'];
CSF_prob_file = [cleaned_dir, '/../region_masks/CSF_prob.nii.gz'];
Edge_prob_file = [cleaned_dir, '/../region_masks/Edge_prob.nii.gz'];
Outbrain_prob_file = [cleaned_dir, '/../region_masks/OutbrainOnly_prob.nii.gz'];
Suscept_prob_file = [cleaned_dir, '/../region_masks/Susceptibility_prob.nii.gz'];
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
[denoised_GM, denoised_GM_mean, ~] = getData_sameregion(denoised_file_data, funcmask, GM_prob_file);
[denoised_GM_NotGM, denoised_NotGM_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, NotGM_prob_file);
[denoised_GM_Edge, denoised_Edge_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Edge_prob_file);
[denoised_GM_Outbrain, denoised_Outbrain_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Outbrain_prob_file);
[denoised_GM_WMCSF, denoised_WMCSF_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, WMCSF_prob_file);
[denoised_GM_CSF, denoised_CSF_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, CSF_prob_file);
[denoised_GM_Suscept, denoised_Suscept_GM, ~, ~, ~, ~] = getData(denoised_file_data, funcmask, GM_prob_file, Suscept_prob_file);

% testing with doing more of the same region
[denoised_NotGM, ~, ~] = getData_sameregion(denoised_file_data, funcmask, NotGM_prob_file);
[denoised_Edge, ~, ~] = getData_sameregion(denoised_file_data, funcmask, Edge_prob_file);
[denoised_Outbrain, ~, ~] = getData_sameregion(denoised_file_data, funcmask, Outbrain_prob_file);
[denoised_WMCSF, ~, ~] = getData_sameregion(denoised_file_data, funcmask, WMCSF_prob_file);
[denoised_CSF, ~, ~] = getData_sameregion(denoised_file_data, funcmask, CSF_prob_file);
[denoised_Suscept, ~, ~] = getData_sameregion(denoised_file_data, funcmask, Suscept_prob_file);

% Now compute relevant correlations, and some others that we don't use, but
% one could make use of if they desired
fprintf('Calculating Relevant Correlations\n')

% create permutations for regions, subtract a little bit because sometimes
% something funky happens and the sizes do not fully align with Edge
% Voxels. Very close (probably related to thresholding values and minor 
% smoothing meaning some voxels that are all 0s for one, are not for the other.)

% For GM compared to other region
GM_randperm = randperm(size(denoised_GM, 1));
GM_Edge_randperm = randperm(size(denoised_GM_Edge, 1));
GM_Outbrain_randperm = randperm(size(denoised_GM_Outbrain, 1));
GM_WMCSF_randperm = randperm(size(denoised_GM_WMCSF, 1));
GM_CSF_randperm = randperm(size(denoised_GM_CSF, 1));
GM_Suscept_randperm = randperm(size(denoised_GM_Suscept, 1));
GM_NotGM_randperm = randperm(size(denoised_GM_NotGM, 1));

% For other region compared to GM
NotGM_GM_randperm = randperm(size(denoised_NotGM_GM, 1));
Outbrain_GM_randperm = randperm(size(denoised_Outbrain_GM, 1));
Edge_GM_randperm = randperm(size(denoised_Edge_GM, 1));
WMCSF_GM_randperm = randperm(size(denoised_WMCSF_GM, 1));
CSF_GM_randperm = randperm(size(denoised_CSF_GM, 1));
Suscept_GM_randperm = randperm(size(denoised_Suscept_GM, 1));

% for region compared to itself:
NotGM_randperm = randperm(size(denoised_NotGM, 1));
Outbrain_randperm = randperm(size(denoised_Outbrain, 1));
Edge_randperm = randperm(size(denoised_Edge, 1));
WMCSF_randperm = randperm(size(denoised_WMCSF, 1));
CSF_randperm = randperm(size(denoised_CSF, 1));
Suscept_randperm = randperm(size(denoised_Suscept, 1));

% Now get all the relevant correlations:
perms = 200;
Edge_GM_corr = createHistData(denoised_GM_Edge, denoised_Edge_GM, GM_Edge_randperm, Edge_GM_randperm, perms);
FD_GM_corr = createHistConf1DData(denoised_GM, confounds_fd, GM_randperm, 19900); % make it same sampling as others
DVARS_GM_corr = createHistConf1DData(denoised_GM, confounds_dvars, GM_randperm, 19900);
Outbrain_GM_corr = createHistData(denoised_GM_Outbrain, denoised_Outbrain_GM, GM_Outbrain_randperm, Outbrain_GM_randperm, perms);
WMCSF_GM_corr = createHistData(denoised_GM_WMCSF, denoised_WMCSF_GM, GM_WMCSF_randperm, WMCSF_GM_randperm, perms);
CSF_GM_corr = createHistData(denoised_GM_CSF, denoised_CSF_GM, GM_CSF_randperm, CSF_GM_randperm, perms);
NotGM_GM_corr = createHistData(denoised_GM_NotGM, denoised_NotGM_GM, GM_NotGM_randperm, NotGM_GM_randperm, perms);
GM_GM_corr = createHistData(denoised_GM, denoised_GM, GM_randperm, GM_randperm, perms);
Suscept_GM_corr = createHistData(denoised_GM_Suscept, denoised_Suscept_GM, GM_Suscept_randperm, Suscept_GM_randperm, perms);

% Getting more of the same regions
Edge_Edge_corr = createHistData(denoised_Edge, denoised_Edge, Edge_randperm, Edge_randperm, perms);
Outbrain_Outbrain_corr = createHistData(denoised_Outbrain, denoised_Outbrain, Outbrain_randperm, Outbrain_randperm, perms);
WMCSF_WMCSF_corr = createHistData(denoised_WMCSF, denoised_WMCSF, WMCSF_randperm, WMCSF_randperm, perms);
CSF_CSF_corr = createHistData(denoised_CSF, denoised_CSF, CSF_randperm, CSF_randperm, perms);
NotGM_NotGM_corr = createHistData(denoised_NotGM, denoised_NotGM, NotGM_randperm, NotGM_randperm, perms);
Suscept_Suscept_corr = createHistData(denoised_Suscept, denoised_Suscept, Suscept_randperm, Suscept_randperm, perms);

GM_mean = denoised_GM_mean;

end

%%
function [GM, region, Global_Signal, GM_Signal_mean, GM_var, region_var] = getData(funcfile_data, funcmask, GM_prob_file, region_prob_file)
% read in the niftis
region_prob_data = niftiread(region_prob_file);
gm_prob_data = niftiread(GM_prob_file);
func_data = funcfile_data;
funcmask_data = niftiread(funcmask);

gm_mask = logical(gm_prob_data > 0.67 & funcmask_data == 1 & region_prob_data == 0 );
region_mask = logical(region_prob_data > 0.67 & funcmask_data == 1 & gm_prob_data < 0.05 );
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

function [region, region_signal_mean, region_var] = getData_sameregion(funcfile_data, funcmask, region_prob_file)
% read in the niftis
region_prob_data = niftiread(region_prob_file);
func_data = funcfile_data;
funcmask_data = niftiread(funcmask);

% Just do 50% probability when getting general region
region_mask = logical(region_prob_data > 0.67 & funcmask_data == 1);
func_mask = logical(funcmask_data == 1);

% Select the timeseries for GM and region and global signal
region = zeros(sum(region_mask(:)), size(func_data, 4));
for j = 1:size(func_data, 4)
    curr_func = func_data(:,:,:,j);
    region(:,j) = curr_func(region_mask);
end

region_var = var(region, 0, 2);
region_signal_mean = mean(region(region(:,1) ~= 0, :));
region_signal_mean = region_signal_mean - mean(region_signal_mean);
end

function compare_orig_corr = createHistData(compare_GM, orig_region, GM_randperm, region_randperm, perms)
% Should center close to 0 (not globally correlated) if it's
% denoised (e.g., denoised lowered noise)
size(compare_GM)
size(GM_randperm)
size(orig_region)
size(region_randperm)
perms
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