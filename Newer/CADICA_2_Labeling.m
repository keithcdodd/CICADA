%% Script to be run after 1_CADICA_MasksandICAs_MNI.sh
% This will label each component as signal or noise
% It is generous -- it follows the "when in doubt, keep as signal."
% Regardless, applying regression of noise ICs after this will greatly
% decrease noise

% This is to follow the 1st bash script.

% Basically, this script labels ICs with reasonably increased gray matter
% overlap and BOLD frequencies as signal and the rest as noise. It attempts
% to label other ICs by noise type, and also outputs "close_ICs" which
% gives a narrowed list of ICs that might want to be double-checked by the
% user.

% Less is more ideas:
% 1. BOLD frequency power > high frequency power
% 2. GM overlap is largest overlap
% Things that can LOOK like GM, but is not: CSF alongside GM boundaries (CSF>GM),
% sinus/veins (Outbrain > Inbrain), motion/edge that correlates with GM (Edge>GM), Transmedullary
% alongside GM borders (WMCSF>GM), Susceptibility (Susc > 0.5 GM)
% 3. At least 50% survives cluster correction
% 4. dvars^2 correlation (data spikes) is not too high (50% of data
% explanation)

clearvars

%%%%%%%%% set up that user may need to be adjust %%%%%%%%%%%%%%%%%%%
workdirext = './'; % this is for first pass
confoundext = './'; % for first pass

%%%%% Everything Else
TR = 2; % in seconds
numvolumes = 300; % how many TRs/samples
cadicafol = '/home/keithdodd/ExampleDataLocal/CADICA_Updated';
% similar syntax as in the 1_CADICA_MasksandICAs
% subjects = {'102' '103' '105' '106' '107' '108' '109' '112' '114' '115' '116' '117' '121' '122' '125' '126' '128' '129' '130' '132' '134' '135' '136' '138' '139' '140' '142' '143' '144' '145' '146' '147' '149' '153' '156' '157' '158' '160' '161' '164' '165' '168' '169' '171' '172' '173' '174' '176' '178' '179' '181' '184' '185' '186' '187'};
subjects = {'102'};
sessions = {'01'};
ROIs = {'GM' 'Edge' 'Transmedullary' 'CSF', 'Susceptibility', 'Outbrain'};
freqs = {'lowfreq', 'BOLDfreq', 'higherfreq'};
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Ultimately, want ICs to lean towards GM and BOLD frequencies, and away from other
% regions and frequencies.

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
       cd(workdirext)

       % (1) Check Spatial Map Overlap with ROIs
       % (1a) Approximate overlap with mean and numvoxels

       % Full volume not cluster corrected
       fullvolume_noclustering_ICmean = dlmread('ROIcalcs/fullvolume_noclustering_ICmean.txt');
       fullvolume_noclustering_ICnumvoxels = dlmread('ROIcalcs/fullvolume_noclustering_ICnumvoxels.txt');
       fullvolume_noclustering_ICsum = fullvolume_noclustering_ICmean .* fullvolume_noclustering_ICnumvoxels;

       % Full volume after cluster correction
       fullvolume_ICmean = dlmread('ROIcalcs/fullvolume_ICmean.txt');
       fullvolume_ICnumvoxels = dlmread('ROIcalcs/fullvolume_ICnumvoxels.txt');
       fullvolume_ICsum = fullvolume_ICmean .* fullvolume_ICnumvoxels;

       % clustering - this will help measure if there is spatial smoothness
       % and actual clusters in the ICs
       clustering_prop = fullvolume_ICsum ./ fullvolume_noclustering_ICsum;

       % smoothness (min cluster size) -> might be useful
       minclustersize = dlmread('clustering/clustersizes.txt');

       % Outbrain
       Outbrain_ICmean = dlmread('ROIcalcs/Outbrain_ICmean.txt');
       Outbrain_ICnumvoxels = dlmread('ROIcalcs/Outbrain_ICnumvoxels.txt');
       Outbrain_ICsum = Outbrain_ICmean .* Outbrain_ICnumvoxels;

       % Inbrain
       Inbrain_ICmean = dlmread('ROIcalcs/Inbrain_ICmean.txt');
       Inbrain_ICnumvoxels = dlmread('ROIcalcs/Inbrain_ICnumvoxels.txt');
       Inbrain_ICsum = Inbrain_ICmean .* Inbrain_ICnumvoxels;

       % Edge
       Edge_ICmean = dlmread('ROIcalcs/Edge_ICmean.txt');
       Edge_ICnumvoxels = dlmread('ROIcalcs/Edge_ICnumvoxels.txt');
       Edge_ICsum = Edge_ICmean .* Edge_ICnumvoxels;

       % Grey Matter
       GM_ICmean = dlmread('ROIcalcs/GM_ICmean.txt');
       GM_ICnumvoxels = dlmread('ROIcalcs/GM_ICnumvoxels.txt');
       GM_ICsum = GM_ICmean .* GM_ICnumvoxels;

       % White Matter CSF Boundary (for subependymal)
       WMCSF_ICmean = dlmread('ROIcalcs/WMCSF_ICmean.txt');
       WMCSF_ICnumvoxels = dlmread('ROIcalcs/WMCSF_ICnumvoxels.txt');
       WMCSF_ICsum = WMCSF_ICmean .* WMCSF_ICnumvoxels;

       % CSF
       CSF_ICmean = dlmread('ROIcalcs/CSF_ICmean.txt');
       CSF_ICnumvoxels = dlmread('ROIcalcs/CSF_ICnumvoxels.txt');
       CSF_ICsum = CSF_ICmean .* CSF_ICnumvoxels;

       % Susceptibility
       Suscept_ICmean = dlmread('ROIcalcs/Suscept_ICmean.txt');
       Suscept_ICnumvoxels = dlmread('ROIcalcs/Suscept_ICnumvoxels.txt');
       Suscept_ICsum = Suscept_ICmean .* Suscept_ICnumvoxels;

       % How many ICs are we examining?
       ICs = 1:length(GM_ICmean);

       % Grab explained variance percent, just in case!
       IC_exp_var = dlmread('ROIcalcs/IC_exp_variance.txt');
       IC_exp_var = IC_exp_var ./ 100; % put it in decimal places

       % (1b) Calculate proportion for each brain region for each IC
       GM_prop = GM_ICsum ./ (Edge_ICsum + GM_ICsum + WMCSF_ICsum + CSF_ICsum + Suscept_ICsum);
       Edge_prop = Edge_ICsum ./ (Edge_ICsum + GM_ICsum + WMCSF_ICsum + CSF_ICsum + Suscept_ICsum);
       WMCSF_prop = WMCSF_ICsum ./ (Edge_ICsum + GM_ICsum + WMCSF_ICsum + CSF_ICsum + Suscept_ICsum);
       CSF_prop = CSF_ICsum ./ (Edge_ICsum + GM_ICsum + WMCSF_ICsum + CSF_ICsum + Suscept_ICsum);
       Suscept_prop = Suscept_ICsum ./ (Edge_ICsum + GM_ICsum + WMCSF_ICsum + CSF_ICsum + Suscept_ICsum);
       Outbrain_prop = (Outbrain_ICsum) ./ (Inbrain_ICsum + Outbrain_ICsum);
       % Outbrain_prop = (Outbrain_ICsum + WMCSF_ICsum) ./ (Inbrain_ICsum + Outbrain_ICsum);

       ROI_props = [GM_prop, Edge_prop, WMCSF_prop, CSF_prop, Suscept_prop, Outbrain_prop];
       ROI_props(isinf(ROI_props)|isnan(ROI_props)) = 0; % replace nonsensicals
       ROI_props_table = array2table(ROI_props, 'VariableNames', ROIs);

       % (2) Check Power Frequency Analysis in low, BOLD (0.008-0.15), and high
       % frequencies
       % (2a) Grab time series and calculate power and proportions in each
       % frequency range of interest
       N = T/dt;
       F = 1/dt;
       df = 1/T;
       N_freq = N/2 + 1;
       lower_phys_cutoff = round(0.008 / df) + 1; % start at freq 0 at position 1
       higher_phys_cutoff = round(0.15 / df) + 1;
       cush = 1; % to avoid overlap, give an index of padding

       % set size of arrays explicitely to save computation time
       ts = zeros(numvolumes, length(ICs));
       f_all_power = zeros(length(ICs), 1);
       BOLDfreqIC = zeros(length(ICs), 1);
       lowfreqIC = zeros(length(ICs), 1);
       highfreqIC = zeros(length(ICs), 1);

       for i=1:length(ICs)
           % grab time spectrum
           file=strcat('melodic/report/t', num2str(i), '.txt');
           ts(:,i)=dlmread(file);

           % calculate power spectrum as evidenced in matlab tutorials
           t = ts(:,i);
           Y = fft(t);
           P2 = abs((Y.^2)/N); % power
           P1 = P2(1:N/2+1);
           f = F*(0:(N/2))/N; % to chart what frequencies we are at

           % capture power range of full spectrum, low, BOLD spectrum, and high Hz
           f_all=trapz(P1);
           f_signal=trapz(P1((lower_phys_cutoff+1):(higher_phys_cutoff-1)));
           f_lowfreq=trapz(P1(1:lower_phys_cutoff));
           f_highfreq=trapz(P1((higher_phys_cutoff+1):length(f)));

           % calculate proportions of each (low, BOLD spectrum, high)
           f_all_power(i) = f_all; % f_all is the same for all though
           BOLDfreqIC(i)=f_signal;
           lowfreqIC(i)=f_lowfreq;
           highfreqIC(i)=f_highfreq;
       end

       % (2b) Calculate actual measured proportions.
       lowfreqIC_prop = lowfreqIC ./ (lowfreqIC+BOLDfreqIC+highfreqIC);
       BOLDfreqIC_prop = BOLDfreqIC ./ (lowfreqIC+BOLDfreqIC+highfreqIC);
       highfreqIC_prop = highfreqIC ./ (lowfreqIC+BOLDfreqIC+highfreqIC);

       freq_props = [lowfreqIC_prop, BOLDfreqIC_prop, highfreqIC_prop];
       freq_props_table = array2table(freq_props, 'VariableNames', freqs);

       % (3) Correlation to data spikes, using dvars
       confound_place = [confoundext, 'confounds_timeseries.csv'];
       allconfounds = readtable(confound_place);
       confounds_9param = table2array(allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'white_matter', 'csf', 'global_signal'}));
       confounds_dvars = table2array(allconfounds(:,{'dvars'}));
       dvars_corr = corr(confounds_dvars(2:end).^4, diff(ts).^4)'.^2;

       % save standard 9Param for comparison:
       writematrix(confounds_9param, [confoundext, '9p_regressors.txt'], 'Delimiter', ' ')

       % (4) Putting it all together!

       % Make a normalized array so it is easier to follow in the future
       feature_array = [ROI_props, freq_props, dvars_corr, clustering_prop];
       factor_names = {'GM', 'Edge', 'Transmed', 'CSF', 'Suscept', 'Outbrain', 'low freq', 'BOLD freq', 'high freq', 'Spike corr', 'Clustering prop'};
       feature_table = array2table(feature_array, 'VariableNames', factor_names);

       % Types of Noise to remove - use a for loop to find ones that are
       % close to being classified as signal or noise
       noise_labelling = zeros(length(ICs), 3);
       d = [-0.05, 0.05, 0];

       for k = 1:3
           curr_GM_prop = GM_prop + d(k);

           % Edge overlap
           Edge_indices = Edge_prop > (GM_prop+d(k));
           Edge_ICs = ICs(Edge_indices);

           % WMCSF instead of GM:
           WMCSF_indices = WMCSF_prop > (GM_prop+d(k));
           WMCSF_ICs = ICs(WMCSF_indices);

           % CSF or internal arteries instead of GM:
           CSF_indices = CSF_prop > (GM_prop+d(k));
           CSF_ICs = ICs(CSF_indices);

           % Susceptibility overlap more than GM
           Suscept_indices = Suscept_prop > (GM_prop+d(k));
           Suscept_ICs = ICs(Suscept_indices);

           % Outbrain instead of inbrain
           Outbrain_indices = Outbrain_prop > (0.5+d(k));
           Outbrain_ICs = ICs(Outbrain_indices);
    
           % Too low of frequency, including susceptibility
           Lowfreq_indices = lowfreqIC_prop > (BOLDfreqIC_prop+d(k));
           Lowfreq_ICs = ICs(Lowfreq_indices);

           % Too high of frequency, including CSF and some motion
           Highfreq_indices = highfreqIC_prop > (BOLDfreqIC_prop+d(k));
           Highfreq_ICs = ICs(Highfreq_indices);

           % Too much spike correlation
           Spikecorr_indices = dvars_corr > (0.5+d(k));
           Spikecorr_ICs = ICs(Spikecorr_indices);

           % Not enough clustering
           Lowcluster_indices = (clustering_prop+d(k)) < 0.5;
           Lowcluster_ICs = ICs(Lowcluster_indices);

           % GM then is anything not labeled as noise
           noise_labelling(:,k) = logical(CSF_indices + Edge_indices + WMCSF_indices + Suscept_indices ...
                + Outbrain_indices + Lowfreq_indices + Highfreq_indices + Spikecorr_indices + Lowcluster_indices);
       end
       
       signal_labelling = ~noise_labelling;

       % Create final noise and signal ICs and indices, if GM_prop>0.5,
       % label as signal regardless
       noise_indices = logical(noise_labelling(:,3) & GM_prop < 0.5);
       noise_ICs = ICs(noise_indices);
       signal_indices = logical(~noise_indices);
       signal_ICs = ICs(signal_indices);
       
       % Create an easy to read table of classification labels
       classification_names = {'IC Number', 'GM-like', 'Edge-like', 'WMCSF-like', 'CSF-like', ...
           'Suscept-like', 'Outbrain-like', 'Low Freq', 'High Freq', 'Spike Correlated', 'Not Clustered'};
       classification_array = [ICs', signal_indices, Edge_indices, WMCSF_indices, ...
           CSF_indices, Suscept_indices, Outbrain_indices, Lowfreq_indices, ...
           Highfreq_indices, Spikecorr_indices, Lowcluster_indices];
       classification_table = array2table(classification_array, 'VariableNames', classification_names);

       % Calculate ones that may be close to being GM or not. We just move
       % the sliding bars a lot, and can examine current GM labels if their
       % GM_prop barely makes the cut. Also if something is labeled both as
       % GM and/or noise.
       close_indices = (noise_labelling(:,1) == 1 & signal_labelling(:,3) == 1) | ...
           (signal_labelling(:,2) == 1 & noise_labelling(:,3) == 1);
       close_ICs = ICs(close_indices)';

       % Make an array and table showing current identification of close ICs
       close_IC_checker = array2table([close_ICs, signal_indices(close_ICs)], 'VariableNames', {'Close IC', 'Labeled as GM Signal?'});
       close_IC_classification = classification_table(close_ICs, :);

       %%%%%%%%%%%%%% Keep editing from here
       % Overall evaluation averages for each type
       means_array = [mean(feature_array(signal_indices,:),1); mean(feature_array(Edge_indices,:),1); ...
           mean(feature_array(WMCSF_indices,:),1); mean(feature_array(CSF_indices,:),1); ...
           mean(feature_array(Suscept_indices,:),1); mean(feature_array(Outbrain_indices,:),1); ...
           mean(feature_array(Lowfreq_indices,:),1); mean(feature_array(Highfreq_indices,:),1); ...
           mean(feature_array(Spikecorr_indices,:),1); mean(feature_array(Lowcluster_indices, :),1);]';
       std_array = [std(feature_array(signal_indices,:),0,1); std(feature_array(Edge_indices,:),0,1); ...
           std(feature_array(WMCSF_indices,:),0,1); std(feature_array(CSF_indices,:),0,1); ...
           std(feature_array(Suscept_indices,:),0,1); std(feature_array(Outbrain_indices,:),0,1); ...
           std(feature_array(Lowfreq_indices,:),0,1); std(feature_array(Highfreq_indices,:),0,1); ...
           std(feature_array(Spikecorr_indices,:),0,1); std(feature_array(Lowcluster_indices, :),0,1);]';
       averages_table = array2table(means_array, 'VariableNames', classification_names(:,2:end), ...
           'RowNames', factor_names);
       std_table = array2table(std_array, 'VariableNames', classification_names(:,2:end), ...
           'RowNames', factor_names);

       % (3b) Make it easy to compare groups
       eval_names = {'IC Number', 'GM', 'Edge', 'Subependymal', 'CSF', 'Suscept', 'Outbrain', 'low freq' 'BOLDfreq', 'HighFreq', 'Spike_Corr', 'Clustering_Prop'};
       eval_table = table(ICs', feature_array(:,1), feature_array(:,2), feature_array(:,3), ...
           feature_array(:,4), feature_array(:,5), feature_array(:,6), feature_array(:,7), ...
           feature_array(:,8), feature_array(:,9), feature_array(:,10), feature_array(:,11), ...
           'VariableNames', eval_names);
       eval_table_close_ICs = eval_table(close_ICs, :);

       % (3d) Compare before and after selection
       before = sum(feature_array .* IC_exp_var);
       IC_exp_var_signal = IC_exp_var(signal_ICs);
       IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
       after = sum(feature_array(signal_ICs,:) .* IC_exp_var_signal);
       compare_cleaning = array2table([before', after'], 'RowNames', eval_names(2:end), 'VariableNames', {'Before', 'After'});

       % (5c) Export variables later use in next steps (e.g., fsl_regfilt)
       writematrix(noise_ICs, 'Noise_dist_ICs.csv')
       writematrix(signal_ICs, 'Signal_dist_ICs.csv')

       % and then make noise components if you want to do aggressive denoising like in CONN (regression)
       mixing_matrix = dlmread('./melodic/melodic_mix');
       noise_dist_covariates = mixing_matrix(:, noise_ICs);
       save('CADICA_Noise_dist.mat', 'noise_dist_covariates')

       % save a matrix of relevant variables if you want to examine later!
       save('DecisionVariables.mat', 'ICs', 'classification_table',...
           'eval_table', 'signal_ICs', 'noise_ICs')
    end
end


% Now, go ahead and compare results and labeling to the FSL outputs in
% Melodic (check report) - you can likely just look at the close ICs.
% You can open the related files too in fsleyes to
% help see the decisions. We suggest overlaying fullvolICA_adj on the
% mnitemplate
