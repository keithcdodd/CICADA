
% set up that may need to be adjusted:
TR = 2;
T = 600;
dt = TR;
curr_wd=pwd;
cd(curr_wd)
addpath(curr_wd);

allconfounds = readtable('confounds_timeseries');
confounds_movement = allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'});
confounds = allconfounds(:,{'global_signal','white_matter','csf', 'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'});


% clusters need to lie mainly in GM, and frequency needs to be mainly in
% BOLD range
GMICmean = dlmread('./GMICmean.txt');
WMICmean = dlmread('./WMICmean.txt');
CSFICmean = dlmread('./CSFICmean.txt');
outbrainICmean = dlmread('./outbrainICmean.txt');

% number of voxels
GMICnumvoxels = dlmread('./GMICnumvoxels.txt');
WMICnumvoxels = dlmread('./WMICnumvoxels.txt');
CSFICnumvoxels = dlmread('./CSFICnumvoxels.txt');
outbrainICnumvoxels = dlmread('./outbrainICnumvoxels.txt');

% approximate sums
GMICsum = GMICmean .* GMICnumvoxels;
WMICsum = WMICmean .* WMICnumvoxels;
CSFICsum = CSFICmean .* CSFICnumvoxels;
outbrainICsum = outbrainICmean .* outbrainICnumvoxels;

% Prepare IC list and reasons for selection
ICs=[1:length(GMICmean)]; % a list of the IC numbers
reasons = cell(length(ICs), 2);

region_significance_names = {'Outbrain', 'GM', 'WM', 'CSF'};
region_significance_table = table(outbrainICsum, GMICsum, WMICsum, CSFICsum, 'VariableNames', region_significance_names);
region_array = [outbrainICsum, GMICsum, WMICsum, CSFICsum];
represented_regions = (region_array == max(region_array, [], 2)) .* logical(sum(region_array, 2)); % what regions are most represented

% Label what region each IC focuses on, if any
outbrain_indices = logical(represented_regions(:,1));
GM_indices = logical(represented_regions(:,2));
WM_indices = logical(represented_regions(:,3));
CSF_indices = logical(represented_regions(:,4));
noregions_indices = ~logical(sum(represented_regions,2));

ICs_outbrain_dominant = ICs(outbrain_indices);
ICs_GM_dominant = ICs(GM_indices);
ICs_WM_dominant = ICs(WM_indices);
ICs_CSF_dominant = ICs(CSF_indices);
% if some did not survive threshold after smoothing - suggests no strong clusters or sources.
ICs_noregions_dominant = ICs(noregions_indices); 


% Now Calculate Based On time & power Spectrum. We want 0.01-0.1 Hz to be
% dominant
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

frequency_significance_names = {'lowfreq', 'BOLDfreq', 'highfreq'};
frequency_significance_table = table(lowfreqICprop, BOLDfreqICprop, highfreqICprop, 'VariableNames', frequency_significance_names);
frequency_array = [lowfreqICprop, BOLDfreqICprop, highfreqICprop];
represented_frequencies = frequency_array == max(frequency_array, [], 2); % What frequency range is most represented?

% test for frequency dominance - we want BOLD frequency dominance
lowfreq_indices = logical(represented_frequencies(:,1));
BOLDfreq_indices = logical(represented_frequencies(:,2)); 
highfreq_indices = logical(represented_frequencies(:,3));

ICs_lowfreq_dominant = ICs(lowfreq_indices);
ICs_BOLDfreq_dominant = ICs(BOLDfreq_indices);
ICs_highfreq_dominant = ICs(highfreq_indices);


% Look for ICAs that highly correlate to movement
confound_corrs = corr(detrend(ts), [detrend(table2array(confounds))]); % detrend because this moves the mean and impacts true correlation
correlated_significant_indices = logical(sum(abs(confound_corrs) > sqrt(0.5), 2));
ICs_correlated_significant = ICs(correlated_significant_indices);


% Finally, look for and obvious spikey data in the time series. The std of
% the absolute values of the derivative should not be more than the mean of
% it. This is conservative. We don't include the first few timepoints in
% case it is simply settling down, but the rest of the signal is decent. We
% control for that later.
spikey_indices = logical((std(abs(diff(ts(3:end,:)))) ./ mean(abs(diff(ts(3:end,:)))))' > 1);
ICs_spikey = ICs(spikey_indices);


% Create table of reasons selected as noise
noisereasons_names = {'Spikey Timeseries', 'Outbrain', 'WM', 'CSF', 'No Region Dominance', 'Low Frequency', 'High Frequency', 'Correlation to Confounds'};
noisefactors_table = table(spikey_indices, outbrain_indices, WM_indices, ...
    CSF_indices, noregions_indices, lowfreq_indices, ...
    highfreq_indices, correlated_significant_indices, 'VariableNames', noisereasons_names);

% GM dominant and BOLD dominant minus any spikey timeseries of highly
% correlated to confound time series
signal_indices = logical(GM_indices .* BOLDfreq_indices .* ~spikey_indices .* ~correlated_significant_indices);
noise_indices = ~signal_indices;

Signal_ICs = ICs(signal_indices);
Noise_ICs = ICs(noise_indices);
% fsl_regfilt will take in the calculated Noise ICs

% now just need to save noise component indices as .csv! then you can run
% fslregfilt!
writematrix(Noise_ICs, 'CalculatedNoiseICs.csv') % good for nonaggressive denoising
writematrix(Signal_ICs, 'CalculatedSignalICs.csv') % good to check it grabbed ones you are happy with!

% and then lets make noise components so we can do aggressive (regression)
% with them if we want
mixing_matrix = dlmread('./melodic/melodic_mix');
noise_covariates = mixing_matrix(:, Noise_ICs);
save('CADICA_Noise.mat', 'noise_covariates')

% save a matrix of variables if you want to examine later!
save('DecisionVariables.mat', 'ICs', 'region_significance_table', 'ICs_spikey','frequency_significance_table', 'noisefactors_table', 'ICs_GM_dominant', 'ICs_BOLDfreq_dominant', 'Noise_ICs', 'Signal_ICs')
