%% Script to be run after 1_CADICA_MasksandICAs
% This will label each component as signal or noise
% there are a few options for this which are outlined and commented
% throughout this code

clearvars

%%%%%%%%% set up that user may need to be adjust %%%%%%%%%%%%%%%%%%%
TR = 2; % in seconds
numvolumes = 300; % how many TRs/samples
cadicafol = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA';
% similar syntax as in the 1_CADICA_MasksandICAs
subjects = {'103'};
sessions = {'01' '02'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



T = TR * numvolumes;
dt = TR;
addpath(cadicafol)

for j = 1:length(subjects)
    cd(cadicafol)
    currsubjfol = ['sub-', subjects{j}];
    cd(currsubjfol)
    for k = 1:length(sessions)
       cd(cadicafol)
       cd(currsubjfol)
       currsessfol = ['ses-', sessions{k}];
       cd(currsessfol)
       
       % (1) Check Spatial Map Overlap with GM vs WM, CSF, and Out of Brain
       % (1a) Approximate overlap with mean and numvoxels 
       GMICmean = dlmread('./GMICmean.txt');
       WMICmean = dlmread('./WMICmean.txt');
       CSFICmean = dlmread('./CSFICmean.txt');
       outbrainICmean = dlmread('./outbrainICmean.txt');
       GMICnumvoxels = dlmread('./GMICnumvoxels.txt');
       WMICnumvoxels = dlmread('./WMICnumvoxels.txt');
       CSFICnumvoxels = dlmread('./CSFICnumvoxels.txt');
       outbrainICnumvoxels = dlmread('./outbrainICnumvoxels.txt');
       GMICsum = GMICmean .* GMICnumvoxels;
       WMICsum = WMICmean .* WMICnumvoxels;
       CSFICsum = CSFICmean .* CSFICnumvoxels;
       outbrainICsum = outbrainICmean .* outbrainICnumvoxels;

       % (1b) Prepare IC list & create easily readable table
       ICs=[1:length(GMICmean)]; % a list of the IC numbers
       reasons = cell(length(ICs), 2);
       region_significance_names = {'Outbrain', 'GM', 'WM', 'CSF'};
       region_significance_table = table(outbrainICsum, GMICsum, WMICsum, CSFICsum, 'VariableNames', region_significance_names);

       % (1c) Calculate Potential Signal vs Noise ICs Based on Spatial Map
       region_array = [outbrainICsum, GMICsum, WMICsum, CSFICsum];
       region_prop = region_array ./ sum(region_array, 2);
       region_prop(isnan(region_prop)) = 0; % convert potential NaN to 0

       % Less aggressive (default): Based on what is most representative
       represented_regions = region_prop == max(region_prop, [], 2);
       % More aggressive alternative: region overlap must hit majority.
       % e.g.:
       % represented_regions = region_prop > 0.5


       % (1d) Label what region each IC focuses on, if any
       outbrain_indices = logical(represented_regions(:,1));
       GM_indices = logical(represented_regions(:,2));
       WM_indices = logical(represented_regions(:,3));
       CSF_indices = logical(represented_regions(:,4));
       noregions_indices = ~logical(sum(represented_regions,2));
       
       ICs_outbrain_dominant = ICs(outbrain_indices);
       ICs_GM_dominant = ICs(GM_indices);
       ICs_WM_dominant = ICs(WM_indices);
       ICs_CSF_dominant = ICs(CSF_indices);
       ICs_noregions_dominant = ICs(noregions_indices);

       % (1e) Testing for more aggressive alternative so you can compare
       GM_majority_indices = region_prop(:,2) > 0.5;
       ICs_GM_majority = ICs(GM_majority_indices);
       
       % (2) Check Power Frequency Analysis in low, BOLD (0.01-0.1), and high
       % frequencies
       % (2a) Grab time series and calculate power and proportions in each
       % frequency range of interest
       N = T/dt;
       F = 1/dt;
       df = 1/T;
       N_freq = N/2 + 1;
       lower_phys_cutoff = 0.01 / df + 1; % start at freq 0 at position 1
       higher_phys_cutoff = 0.1 / df + 1;

       % set size of arrays explicitely to save computation time
       ts = zeros(numvolumes, length(ICs));
       f_all_power = zeros(length(ICs), 1);
       BOLDfreqICprop = zeros(length(ICs), 1);
       lowfreqICprop = zeros(length(ICs), 1);
       highfreqICprop = zeros(length(ICs), 1);
       for i=1:length(ICs)   
           % grab time spectrum
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
           f_all_power(i) = f_all;
           BOLDfreqICprop(i)=f_signal/f_all;
           lowfreqICprop(i)=f_lowfreq/f_all;
           highfreqICprop(i)=f_highfreq/f_all;
        end
        
       % (2b) Create easily readable table
       frequency_significance_names = {'lowfreq', 'BOLDfreq', 'highfreq'};
       frequency_significance_table = table(lowfreqICprop, BOLDfreqICprop, highfreqICprop, 'VariableNames', frequency_significance_names);
        
       % (2c) Calculate Potential Signal vs Noise ICs Based on Power
       % Frequency
       frequency_array = [lowfreqICprop, BOLDfreqICprop, highfreqICprop];
       frequency_prop = frequency_array ./ sum(frequency_array, 2);

       % Less aggressive (default): Based on what is most representative 
       represented_frequencies = frequency_prop == max(frequency_prop, [], 2);
       % More aggressive alternative: BOLD power must hit majority. e.g.
       % represented_frequencies = frequency_prop > 0.5

       % (2d) Label what frequency range each IC focuses on, if any 
       % test for frequency dominance - we want BOLD frequency dominance
       lowfreq_indices = logical(represented_frequencies(:,1));
       BOLDfreq_indices = logical(represented_frequencies(:,2)); 
       highfreq_indices = logical(represented_frequencies(:,3));
        
       ICs_lowfreq_dominant = ICs(lowfreq_indices);
       ICs_BOLDfreq_dominant = ICs(BOLDfreq_indices);
       ICs_highfreq_dominant = ICs(highfreq_indices);

       % (2e) Testing for more aggressive alternative so you can compare
       BOLD_majority_indices = frequency_prop(:,2) > 0.5;
       ICs_BOLD_majority = ICs(BOLD_majority_indices);
       
       % (2f) Calculate mean proportion between (1) and (2) in case you
       % want to do medium aggressive selection later
       GM_BOLD_prop = (region_prop(:,2) + frequency_prop(:,2)) ./ 2;
       % If mean propotion is greater than 50%, then more certain that it
       % is signal
       GM_BOLD_majority_indices = logical(GM_BOLD_prop > 0.5);
       ICs_GM_BOLD_majority = ICs(GM_BOLD_majority_indices);

       % (3) Identify Highly Spikey time series data
       %  If the std is more than the mean of the absolute values, then it is very spikey data (conservative approach). 
       % We don't include the first few timepoints in case it is simply
       % settling down - either delete first couple volumes or control for
       % first few volumes in a GLM - "control for task effect"
       spikey_indices = logical((std(abs(diff(ts(3:end,:)))) ./ mean(abs(diff(ts(3:end,:)))))' > 1);
       ICs_spikey = ICs(spikey_indices);

       % (4) Identify ICs highly correlated to common confounds
       % (4a) Grab confound information
       allconfounds = readtable('confounds_timeseries');
       confounds_general = allconfounds(:,{'global_signal','white_matter','csf', 'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'});
       
       % (4b) Look for ICAs that highly correlate to common confounds
       % detrend because of otherwise moving means which can hide a strong correlation
       confound_corrs = corr(detrend(ts), [detrend(table2array(confounds_general))]); 
       
       % (4c) Only label ICs where over 50% of the variance can be accounted for
       % due to a confound (very conservative)
       correlated_significant_indices = logical(sum(abs(confound_corrs) > sqrt(0.5), 2));
       ICs_correlated_significant = ICs(correlated_significant_indices);
       
       % (5) Putting it all together!
       % (5a) Create an easy to read table of noise reasons
       noisereasons_names = {'Spikey Timeseries', 'Outbrain', 'WM', 'CSF', 'No Region Dominance', 'Low Frequency', 'High Frequency', 'Correlation to Confounds'};
       noisefactors_table = table(spikey_indices, outbrain_indices, WM_indices, ...
           CSF_indices, noregions_indices, lowfreq_indices, ...
           highfreq_indices, correlated_significant_indices, 'VariableNames', noisereasons_names);

       % (5b) Label Signal ICs vs Noise ICs
       % Less Selective (default): (Labels less signal as noise)
       signal_lowsel_indices = logical(GM_indices .* BOLDfreq_indices .* ~spikey_indices .* ~correlated_significant_indices);
       % Medium Selective for comparison: incorporates (2f)
       signal_medsel_indices = logical(GM_indices .* BOLDfreq_indices .* GM_BOLD_majority_indices .* ~spikey_indices .* ~correlated_significant_indices);
       % More Selective for comparison: (Labels less noise as signal)
       signal_highsel_indices = logical(GM_majority_indices .* BOLD_majority_indices .* ~spikey_indices .* ~correlated_significant_indices);

       Signal_lowsel_ICs = ICs(signal_lowsel_indices);
       Signal_medsel_ICs = ICs(signal_medsel_indices);
       Signal_highsel_ICs = ICs(signal_highsel_indices);
       Noise_lowsel_ICs = ICs(~signal_lowsel_indices);
       Noise_medsel_ICs = ICs(~signal_medsel_indices);
       Noise_highsel_ICs = ICs(~signal_highsel_indices);

       % there may be times in difficult data where higher selection gives
       % 0 ICs as signal. We want to be able to handle that instead of
       % getting an error in write matrix. So set those to just 0.
       if sum(signal_highsel_indices) < 1
           Signal_highsel_ICs = 0;
           if sum(signal_medsel_indices) < 1
               Signal_medsel_ICs = 0;
           end
       end

       % (5c) Export variables later use in next steps (e.g., fsl_regfilt)
       writematrix(Noise_lowsel_ICs, 'Noise_lowsel_ICs.csv')
       writematrix(Noise_medsel_ICs, 'Noise_medsel_ICs.csv')
       writematrix(Noise_highsel_ICs, 'Noise_highsel_ICs.csv')
       writematrix(Signal_lowsel_ICs, 'Signal_lowsel_ICs.csv')
       writematrix(Signal_medsel_ICs, 'Signal_medsel_ICs.csv')
       writematrix(Signal_highsel_ICs, 'Signal_highsel_ICs.csv')
        
       % and then make noise components if you want to do aggressive denoising (regression)
       mixing_matrix = dlmread('./melodic/melodic_mix');
       noise_lowsel_covariates = mixing_matrix(:, Noise_lowsel_ICs);
       noise_medsel_covariates = mixing_matrix(:, Noise_medsel_ICs);
       noise_highsel_covariates = mixing_matrix(:, Noise_highsel_ICs);
       save('CADICA_Noise_lowsel.mat', 'noise_lowsel_covariates')
       save('CADICA_Noise_medsel.mat', 'noise_medsel_covariates')
       save('CADICA_Noise_highsel.mat', 'noise_highsel_covariates')

       % save a matrix of relevant variables if you want to examine later!
       save('DecisionVariables.mat', 'ICs', 'region_significance_table', 'ICs_spikey','frequency_significance_table', ...
           'noisefactors_table', 'ICs_GM_dominant', 'ICs_BOLDfreq_dominant', 'ICs_GM_BOLD_majority', ...
           'ICs_GM_majority', 'ICs_BOLD_majority', 'ICs_correlated_significant','Signal_lowsel_ICs', 'Signal_medsel_ICs', 'Signal_highsel_ICs')
    end
end


% Now, go ahead and compare results and labeling to the FSL outputs in
% Melodic (check report). You can open the related files too in fsleyes to
% help see the decisions. We suggest overlaying fullvolICA_adj on the
% mnitemplate


