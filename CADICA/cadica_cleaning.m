
clearvars

% set up that may need to be adjusted:
TR = 2;
T = 600;
dt = TR;
curr_wd=pwd;
cd(curr_wd)
addpath(curr_wd);
% get relevant folders for session & subject
out = regexp(pwd,'/','split');
session = out{end};
subject = out{end-1};


% load confound timeseries file and grab what you need/want trans_x, trans_y, trans_z, rot_x,
% rot_y, rot_z, csf, white_matter, global_signal, take the normalized above
% 1.65 for each and add those on as nosiest and noisy signals
allconfounds = readtable('confounds_timeseries');
confounds = allconfounds(:,{'global_signal','white_matter','csf', 'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'csf_wm'});
confounds_movement = allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'});
translations = confounds(:,{'trans_x', 'trans_y', 'trans_z'});
rotations = confounds(:,{'rot_x', 'rot_y', 'rot_z'});
abs_trans = sqrt(sum(table2array(translations) .^2, 2));
abs_rotations = sqrt(sum(table2array(rotations) .^2, 2));
framewise_displacement = allconfounds(2:end, {'framewise_displacement'});


% start with weights for GM, WM, CSF, inbrain, outbrain
% gather means
GMprobmean = dlmread('./GMprobmean.txt');
WMprobmean = dlmread('./WMprobmean.txt');
CSFprobmean = dlmread('./CSFprobmean.txt');
inbrainprobmean = dlmread('./inbrainprobmean.txt');
outbrainprobmean = dlmread('./outbrainprobmean.txt');
% gather numvoxels
GMprobnumvoxels = dlmread('./GMprobnumvoxels.txt');
WMprobnumvoxels = dlmread('./WMprobnumvoxels.txt');
CSFprobnumvoxels = dlmread('./CSFprobnumvoxels.txt');
inbrainprobnumvoxels = dlmread('./inbrainprobnumvoxels.txt');
outbrainprobnumvoxels = dlmread('./outbrainprobnumvoxels.txt');
% calculate approximate sums/weights
GMprobsum = GMprobmean .* GMprobnumvoxels;
WMprobsum = WMprobmean .* WMprobnumvoxels;
CSFprobsum = CSFprobmean .* CSFprobnumvoxels;
inbrainprobsum = inbrainprobmean .* inbrainprobnumvoxels;
outbrainprobsum = outbrainprobmean .* outbrainprobnumvoxels;
outGMprobsum = WMprobsum + CSFprobsum + outbrainprobsum;

% and now calculated expected proportions
inbrainprop = inbrainprobsum / (GMprobsum + WMprobsum + CSFprobsum + outbrainprobsum);
outbrainprop = outbrainprobsum / (GMprobsum + WMprobsum + CSFprobsum + outbrainprobsum);
GMprop = GMprobsum / (GMprobsum + WMprobsum + CSFprobsum + outbrainprobsum);
WMprop = WMprobsum / (GMprobsum + WMprobsum + CSFprobsum + outbrainprobsum);
CSFprop = CSFprobsum / (GMprobsum + WMprobsum + CSFprobsum + outbrainprobsum);
outGMprop = outGMprobsum / (GMprobsum + WMprobsum + CSFprobsum + outbrainprobsum); 

% now move onto actual sampled results
% means
GMICmean = dlmread('./GMICmean.txt');
WMICmean = dlmread('./WMICmean.txt');
CSFICmean = dlmread('./CSFICmean.txt');
inbrainICmean = dlmread('./inbrainICmean.txt');
outbrainICmean = dlmread('./outbrainICmean.txt');
% number of voxels
GMICnumvoxels = dlmread('./GMICnumvoxels.txt');
WMICnumvoxels = dlmread('./WMICnumvoxels.txt');
CSFICnumvoxels = dlmread('./CSFICnumvoxels.txt');
inbrainICnumvoxels = dlmread('./inbrainICnumvoxels.txt');
outbrainICnumvoxels = dlmread('./outbrainICnumvoxels.txt');
% approximate sums
GMICsum = GMICmean .* GMICnumvoxels;
WMICsum = WMICmean .* WMICnumvoxels;
CSFICsum = CSFICmean .* CSFICnumvoxels;
inbrainICsum = inbrainICmean .* inbrainICnumvoxels;
outbrainICsum = outbrainICmean .* outbrainICnumvoxels;
outGMICsum = WMICsum + CSFICsum + outbrainICsum;

% calculate sampled proportions
inbrainICprop = inbrainICsum ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum);
outbrainICprop = outbrainICsum ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum);
GMICprop = GMICsum ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum);
WMICprop = WMICsum ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum);
CSFICprop = CSFICsum ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum);
outGMICprop = outGMICsum ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum);

% Prepare IC list and reasons for selection
ICs=[1:length(inbrainICmean)]; % a list of the IC numbers
reasons = cell(length(ICs), 2);

% put the relevant proportions all together for later reference
% prop_array = [outeredgeIC_prop, outsidebrainIC_prop, insidebrainIC_prop, GMIC_prop, WMIC_prop, CSFIC_prop];
propICarray = [outbrainICprop, inbrainICprop, outGMICprop, GMICprop, WMICprop, CSFICprop];
proparray = [outbrainprop, inbrainprop, outGMprop, GMprop, WMprop, CSFprop];

% test for proportional changes of outbrain, inbrain, GM, WM, and CSF
% first in vs outbrain
GM_z_prop = (GMICprop - GMprop) ./ sqrt((GMprop.*(1-GMprop)) ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum));
WM_z_prop = (WMICprop - WMprop) ./ sqrt((WMprop.*(1-WMprop)) ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum));
CSF_z_prop = (CSFICprop - CSFprop) ./ sqrt((CSFprop.*(1-CSFprop)) ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum));
outbrain_z_prop = (outbrainICprop - outbrainprop) ./ sqrt((outbrainprop.*(1-outbrainprop)) ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum));
outGM_z_prop = (outGMICprop - outGMprop) ./ sqrt((outGMprop.*(1-outGMprop)) ./ (GMICsum + WMICsum + CSFICsum + outbrainICsum));

z_prop_regions = [outbrain_z_prop, outGM_z_prop, GM_z_prop, WM_z_prop, CSF_z_prop];

region_significance_names = {'Outbrain', 'OutGM', 'GM', 'WM', 'CSF'};
region_significance_table = table(outbrain_z_prop, outGM_z_prop, GM_z_prop, WM_z_prop, CSF_z_prop, 'VariableNames', region_significance_names);

% test where GM relative proportion is significant
GM_significant_indices = normalize(GM_z_prop) > 1.65;
ICs_GM_significant = ICs(GM_significant_indices);

% test for more significant WM signal within the brain
WM_significant_indices = normalize(WM_z_prop) > 1.65;
ICs_WM_significant = ICs(WM_significant_indices);

% test for more significant CSF signal within the brain
CSF_significant_indices = normalize(CSF_z_prop) > 1.65;
ICs_CSF_significant = ICs(CSF_significant_indices);

% test for more significant signal outside the brain
outbrain_significant_indices = normalize(outbrain_z_prop) > 1.65;
ICs_outbrain_significant = ICs(outbrain_significant_indices);

% test for more significant signal outside of GM
outGM_significant_indices = normalize(outGM_z_prop) > 1.65;
ICs_outGM_significant = ICs(outGM_significant_indices);

% calculate indices that did not meet any criteria
no_significantregions_indices = ~(outbrain_significant_indices + outGM_significant_indices + GM_significant_indices + WM_significant_indices + CSF_significant_indices);
ICs_no_significantregions = ICs(no_significantregions_indices);

region_noise_indices = logical(outbrain_significant_indices + WM_significant_indices + CSF_significant_indices);
region_signal_indices = logical(GM_significant_indices .* ~region_noise_indices);
ICs_region_significant = ICs(region_signal_indices);


% Now Calculate Based On Power Spectrum. We want 0.01-0.1 Hz
N = T/dt;
F = 1/dt;
df = 1/T;
N_freq = N/2 + 1;
lower_phys_cutoff = 0.01 / df + 1; % start at freq 0 at position 1
higher_phys_cutoff = 0.1 / df + 1;
for i=1:length(ICs)   
    % time spectrum
    file=strcat('./melodic/report/t', num2str(i), '.txt');
    ts(:,i)=dlmread(file);

    % calculate power spectrum as evidenced in matlab tutorials
    t = ts(:,i);
    Y = fft(t);
    P2 = abs((Y.^2)/N); % power
    P1 = P2(1:N/2+1);
    f = F*(0:(N/2))/N; % to chart what frequencies we are at

    % capture power range of full spectrum, low, BOLD spectrum, and high Hz
    f_all=trapz(P1);
    f_signal=trapz(P1(lower_phys_cutoff:higher_phys_cutoff));
    f_lowfreq=trapz(P1(1:lower_phys_cutoff));
    f_highfreq=trapz(P1(higher_phys_cutoff:length(f)));

    % calculate proportions of each (low, BOLD spectrum, high)
    f_all_power(i,1) = f_all;
    BOLDfreqICprop(i,1)=f_signal/f_all;
    lowfreqICprop(i,1)=f_lowfreq/f_all;
    highfreqICprop(i,1)=f_highfreq/f_all;
end
% save a table for later reference:
fICprops = [lowfreqICprop, BOLDfreqICprop, highfreqICprop];

% Calculate expected proportions for low, BOLD, and high
lowfreqprop = (lower_phys_cutoff - 1) / N_freq;
BOLDfreqprop = (higher_phys_cutoff - lower_phys_cutoff) / N_freq;
highfreqprop = (N_freq - higher_phys_cutoff) / N_freq;
fprops = [lowfreqprop, BOLDfreqprop, highfreqprop];

% test for proportional changes of power for each of the three ranges
lowfreq_z_prop = (lowfreqICprop - lowfreqprop) ./ sqrt((lowfreqprop.*(1-lowfreqprop)) ./ N_freq);
BOLDfreq_z_prop = (BOLDfreqICprop - BOLDfreqprop) ./ sqrt((BOLDfreqprop.*(1-BOLDfreqprop)) ./ N_freq);
highfreq_z_prop = (highfreqICprop - highfreqprop) ./ sqrt((highfreqprop.*(1-highfreqprop)) ./ N_freq);

z_prop_frequencies = [lowfreq_z_prop, BOLDfreq_z_prop, highfreq_z_prop];

% test for when there is significant and more proportional signal increase
% in low frequency power
lowfreq_significant_indices = normalize(lowfreq_z_prop) > 1.65;
ICs_lowfreq_significant = ICs(lowfreq_significant_indices);

% test for when there is significant and more proportional signal increase
% in BOLD frequency power
% sometimes true signal hides well in low frequency power. We high pass
BOLDfreq_significant_indices = normalize(BOLDfreq_z_prop) > 1.65;
ICs_BOLDfreq_significant = ICs(BOLDfreq_significant_indices);

% test for when there is significant and more proportional signal increase
% in high frequency power
highfreq_significant_indices = normalize(highfreq_z_prop) > 1.65;
ICs_highfreq_significant = ICs(highfreq_significant_indices);

% calculate indices that did not meet any criteria
nofreq_significant_indices = ~(lowfreq_significant_indices + highfreq_significant_indices);
ICs_nofreq_significant = ICs(nofreq_significant_indices);

frequency_noise_indices = logical(lowfreq_significant_indices + highfreq_significant_indices);
frequency_signal_indices = logical(BOLDfreq_significant_indices .* ~frequency_noise_indices);
ICs_frequency_significant = ICs(frequency_signal_indices);

% Do overlap of GM and BOLD here. Because good ICs are ones that show GM and BOLD freq together!
GM_BOLD_overlap_indices = normalize(normalize(GM_z_prop) + normalize(BOLDfreq_z_prop)) > 1.65;
ICs_GM_BOLD_overlap = ICs(GM_BOLD_overlap_indices);

% (3) now lets do the correlation stuff with the confounds of interest - remove
% high positive and low negative correlation. 0.7 because thats about 50%
% variance explained by confounding factor
% fd_corrs = corr(ts(1:end-1, :), table2array(framewise_displacement)); % do this separately because of 1 less input
confound_corrs = corr(ts, [table2array(confounds_movement)]);
correlated_significant_indices = logical(sum(abs(confound_corrs) > sqrt(0.5), 2));
ICs_correlated_significant = ICs(correlated_significant_indices);

% (4) Let's also look for time series with spikes or massive drifts in it in general 
spikey_significant_indices = logical(normalize(max(abs(diff(ts))) ./ mean(abs(diff(ts)))) > 1.65)';
ICs_spikey_significant = ICs(spikey_significant_indices);

% (5) And also look for time series with spikes that correlate to motion
motionspike_significant_indices = logical((normalize(sum(abs(diff(abs_trans)) .* abs(diff(ts)))) > 1.65) + ...
    (normalize(sum(abs(diff(abs_rotations)) .* abs(diff(ts)))) > 1.65) + ...
    (normalize(sum(abs(table2array(framewise_displacement)) .* abs(diff(ts)))) > 1.65))';
ICs_motionspike_significant = ICs(motionspike_significant_indices);


% Create table of reasons selected as noise
noisereasons_names = {'Correlated Confounds', 'Spikey Timeseries', 'Motion Spikes', 'Outbrain', 'WM', 'CSF', 'No Region Dominance', 'Low Frequency', 'High Frequency', 'No Frequency Dominance'};
noisefactors_table = table(correlated_significant_indices, spikey_significant_indices, motionspike_significant_indices, outbrain_significant_indices, WM_significant_indices, ...
    CSF_significant_indices, no_significantregions_indices, lowfreq_significant_indices, ...
    highfreq_significant_indices, nofreq_significant_indices, 'VariableNames', noisereasons_names);

% signal will be the multiplication of the 2 indices
% strict_signal_indices = logical(GM_significant_indices .* BOLDfreq_significant_indices .* GM_BOLD_overlap_indices);
noisiest_indices = logical(region_noise_indices + frequency_noise_indices + ...
    correlated_significant_indices + spikey_significant_indices + motionspike_significant_indices);
signal_indices = logical(GM_BOLD_overlap_indices .* ~noisiest_indices); % in case there is any overlap
noisy_indices = ~signal_indices;

Signal_ICs = ICs(signal_indices);
Noisiest_ICs = ICs(noisiest_indices); % if I do noisiest only, then I like aggressive denoising for fsl_regfilt
Noisy_ICs = ICs(noisy_indices); % this might be good for non-aggressive
% fsl_regfilt will take in the calculated Noise ICs

% now just need to save noise component indices as .csv! then you can run
% fslregfilt!
writematrix(Noisy_ICs, 'CalculatedNoisyICs.csv') % good for nonaggressive denoising
writematrix(Noisiest_ICs, 'CalculatedNoisiestICs.csv') % good for aggressive denoising
writematrix(ICs_correlated_significant, 'CalculatedCorrelatedICs.csv') % good for aggressive denoising - will be less than Noisiest
writematrix(Signal_ICs, 'CalculatedSignalICs.csv')

% and then lets make noise components so we can do aggressive (regression)
% with them if we want
mixing_matrix = dlmread('./melodic/melodic_mix');
noisiest_covariates = mixing_matrix(:, Noisiest_ICs);
noisy_covariates = mixing_matrix(:, Noisy_ICs);
save('CADICA_noisiest.mat', 'noisiest_covariates')
save('CADICA_noisy.mat', 'noisy_covariates')

% save a matrix of variables if you want to examine later!
save('DecisionVariables.mat', 'ICs', 'propICarray', 'proparray', 'fICprops', 'fprops', 'z_prop_regions', 'z_prop_frequencies', 'noisefactors_table', 'ICs_correlated_significant', 'ICs_frequency_significant', 'ICs_GM_BOLD_overlap', 'ICs_region_significant', 'Signal_ICs')

