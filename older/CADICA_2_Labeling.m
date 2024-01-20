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

clearvars

%%%%%%%%% set up that user may need to be adjust %%%%%%%%%%%%%%%%%%%

%%%%% First, need to change a couple things depending on 1st or 2nd pass:
workdirext = './'; % this is for first pass
confoundext = './'; % for first pass
%workdirext = 'cleaned'; % this is for second pass
%confoundext = '../'; % for second pass

%%%%% Everything Else:signa
TR = 2; % in seconds
numvolumes = 300; % how many TRs/samples
cadicafol = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA_Updated';
% similar syntax as in the 1_CADICA_MasksandICAs
% subjects = {'102' '103' '105' '106' '107' '108' '109' '112' '114' '115' '116' '117' '121' '122' '125' '126' '128' '129' '130' '132' '134' '135' '136' '138' '139' '140' '142' '143' '144' '145' '146' '147' '149' '153' '156' '157' '158' '160' '161' '164' '165' '168' '169' '171' '172' '173' '174' '176' '178' '179' '181' '184' '185' '186' '187'};
subjects = {'102'};
sessions = {'01'};
ROIs = {'GM' 'Edge' 'Transmedullary' 'CSF', 'Susceptibility'};
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

       % smoothness (min cluster size)
       minclustersize = dlmread('clustering/clustersizes.txt');

       % full volume
       fullvolume_ICmean = dlmread('ROIcalcs/fullvolume_ICmean.txt');
       fullvolume_ICnumvoxels = dlmread('ROIcalcs/fullvolume_ICnumvoxels.txt');
       fullvolume_ICsum = fullvolume_ICmean .* fullvolume_ICnumvoxels;

       % outbrain
       Outbrain_ICmean = dlmread('ROIcalcs/Outbrain_ICmean.txt');
       Outbrain_ICnumvoxels = dlmread('ROIcalcs/Outbrain_ICnumvoxels.txt');
       Outbrain_ICsum = Outbrain_ICmean .* Outbrain_ICnumvoxels;

       % Edge
       Edge_ICmean = dlmread('ROIcalcs/Edge_ICmean.txt');
       Edge_ICnumvoxels = dlmread('ROIcalcs/Edge_ICnumvoxels.txt');
       Edge_ICsum = Edge_ICmean .* Edge_ICnumvoxels;

       % Grey Matter
       GM_ICmean = dlmread('ROIcalcs/GM_ICmean.txt');
       GM_ICnumvoxels = dlmread('ROIcalcs/GM_ICnumvoxels.txt');
       GM_ICsum = GM_ICmean .* GM_ICnumvoxels;

       % White Matter
       WM_ICmean = dlmread('ROIcalcs/WM_ICmean.txt');
       WM_ICnumvoxels = dlmread('ROIcalcs/WM_ICnumvoxels.txt');
       WM_ICsum = WM_ICmean .* WM_ICnumvoxels;

       % White Matter CSF Boundary (for subependymal)
       WMCSF_ICmean = dlmread('ROIcalcs/WMCSF_ICmean.txt');
       WMCSF_ICnumvoxels = dlmread('ROIcalcs/WMCSF_ICnumvoxels.txt');
       WMCSF_ICsum = WMCSF_ICmean .* WMCSF_ICnumvoxels;

       % CSF
       CSF_ICmean = dlmread('ROIcalcs/CSF_ICmean.txt');
       CSF_ICnumvoxels = dlmread('ROIcalcs/CSF_ICnumvoxels.txt');
       CSF_ICsum = CSF_ICmean .* CSF_ICnumvoxels;

       % Inbrain
       Inbrain_ICmean = dlmread('ROIcalcs/Inbrain_ICmean.txt');
       Inbrain_ICnumvoxels = dlmread('ROIcalcs/Inbrain_ICnumvoxels.txt');
       Inbrain_ICsum = Inbrain_ICmean .* Inbrain_ICnumvoxels;

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
       Inbrain_prop = Inbrain_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       Outbrain_prop = Outbrain_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       Edge_prop = Edge_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       GM_prop = GM_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       WM_prop = WM_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       WMCSF_prop = WMCSF_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       CSF_prop = CSF_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);
       Suscept_prop = Suscept_ICsum ./ (Edge_ICsum + GM_ICsum + WM_ICsum + CSF_ICsum);

       ROI_props_nonadj = [GM_prop, Edge_prop, WMCSF_prop, CSF_prop, Suscept_prop]; % we don't need Outbrain or WM
       ROI_props = ROI_props_nonadj ./ sum(ROI_props_nonadj(:,[1:4]),2);
       ROI_props_more = [GM_prop, Outbrain_prop, Edge_prop, WMCSF_prop, CSF_prop, Suscept_prop];
       ROI_props(isinf(ROI_props)|isnan(ROI_props)) = 0; % replace nonsensicals
       ROI_props_table = array2table(ROI_props, 'VariableNames', ROIs);
       Outbrain_prop = ROI_props(:,2) + ROI_props(:,4); % Edge plus CSF

       % Check coverage for ROIs:
       GM_dom_indices = ROI_props(:,1) > Outbrain_prop;
       Edge_dom_indices = (ROI_props(:,2) == max(ROI_props, [], 2) | (ROI_props(:,2) > mean(ROI_props(:,2)) & isoutlier(ROI_props(:,2), "gesd")) );
       WMCSF_dom_indices = (ROI_props(:,3) == max(ROI_props, [], 2) | (ROI_props(:,3) > mean(ROI_props(:,3)) & isoutlier(ROI_props(:,3), "gesd")) );
       CSF_dom_indices = (ROI_props(:,4) == max(ROI_props, [], 2) | (ROI_props(:,4) > mean(ROI_props(:,4)) & isoutlier(ROI_props(:,4), "gesd")) );
       Suscept_dom_indices = (ROI_props(:,5) == max(ROI_props, [], 2) | (ROI_props(:,5) > mean(ROI_props(:,5)) & isoutlier(ROI_props(:,5), "gesd")) );
       Outbrain_dom_indices = ((ROI_props(:,2)+ROI_props(:,4)) > max(ROI_props, [], 2) | ...
           ((ROI_props(:,2)+ROI_props(:,4)) > mean(ROI_props(:,2)+ROI_props(:,4)) & isoutlier(ROI_props(:,2)+ROI_props(:,4), "gesd")) ); % This will be a LOT of them
       

       % (2) Check Power Frequency Analysis in low, BOLD (0.008-0.15), and high
       % frequencies
       % (2a) Grab time series and calculate power and proportions in each
       % frequency range of interest
       N = T/dt;
       F = 1/dt;
       df = 1/T;
       N_freq = N/2 + 1;
       lower_phys_cutoff = round(0.008 / df) + 1; % start at freq 0 at position 1
       higher_phys_cutoff = round(0.1 / df) + 1;
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

       freq_props = [lowfreqIC_prop, BOLDfreqIC_prop, highfreqIC_prop]; %don't need lowfreq
       freq_props_table = array2table(freq_props, 'VariableNames', freqs);

       % (2c) Select freq props based on majority like with GM
       Lowfreq_dom_indices = (freq_props(:,1) == max(freq_props, [], 2) | (freq_props(:,1) > mean(freq_props(:,1)) & isoutlier(freq_props(:,1), "gesd")) );
       BOLDfreq_dom_indices = (freq_props(:,2) == max(freq_props, [], 2) | (freq_props(:,2) > mean(freq_props(:,2)) & isoutlier(freq_props(:,2), "gesd")) );
       Highfreq_dom_indices = (freq_props(:,3) == max(freq_props, [], 2) | (freq_props(:,3) > mean(freq_props(:,3)) & isoutlier(freq_props(:,3), "gesd")) );

       % (3) Correlation to movement spikes
       confound_place = [confoundext, 'confounds_timeseries'];
       allconfounds = readtable(confound_place);
       confounds_motion = allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'framewise_displacement'});
       confounds_motion_array = table2array(confounds_motion);
       confounds_motion_array(1,end) = 0; % because fd starts with NaN
       %confounds_motion_array = normalize(confounds_motion_array).^4;

       % Just run cross correlations, and if any of them explain > 50% of
       % the data, that is bad
       numlags = 2;
       ts_padded = padarray(ts, [numlags 0], 'both');
       corrs = zeros(length(ICs), 7, 4);
       % max_partialcor = zeros(numlags+1, width(table2array(confounds_motion)));
       
       % we step through to do a TR before the motion is recorded, and then
       % 2 TRs after motion is recorded
       for j = 1:(numlags+2)
           for k = 1:width(ts)
            curr_partialcorr = corr(ts_padded(numlags-1+j:(end-numlags-2+j),k), confounds_motion_array);
            corrs(k,1:width(confounds_motion_array),j) = curr_partialcorr;
           end
       end
       
       % So we can find the total summation, for each lag, of timeseries
       % that can be accounted for by motion parameters
       MaxCoeffsLags = max(corrs.^2, [], 3);
       % Then select the max of those lags
       maxCoeffsConfounds = max( MaxCoeffsLags, [], 2);

       motion_corr = maxCoeffsConfounds;

       % (4) Putting it all together!

       % Make a normalized array so it is easier to follow in the future
       feature_array = [ROI_props, freq_props(:,2:3), motion_corr, clustering_prop];
       factor_names = {'GM', 'Edge', 'Transmed', 'CSF', 'Suscept', 'BOLD freq', 'Higher freq', 'Motion Corr', 'Clustering Prop'};
       feature_table = array2table(feature_array, 'VariableNames', factor_names);
       X = normalize(feature_array); % mean is 0, std is 1
       X_table = array2table(X, 'VariableNames', factor_names);

       negative_factors = [ROI_props(:,2:4) ./ ROI_props(:,1), sum(ROI_props(:,2:4),2) ./ ROI_props(:,1)];
       negative_names = {'Edge/GM', 'Transmed/GM', 'CSF/GM', 'OutGM/GM'};
       negative_factor_table = array2table(negative_factors, 'VariableNames', negative_names);

       X_neg = normalize(negative_factors);
       X_neg_table = array2table(X_neg, 'VariableNames',negative_names);
       
       % Bandpassing first helps us distinguish! Not smoothing yet though,
       % retains small ROIs
      

       % Regions outside of GM are not disproportionally represented
       % compared to GM, clustering is at least 50%, and motion corr is not
       % horrendous
       GM_pot_indices = sum(X_neg < 0, 2) == width(X_neg) & clustering_prop > 0.5 & normalize(motion_corr) < 1;
       GM_pot_ICs = ICs(GM_pot_indices);

       % Motion or sudden signal artifact, or Sinus: Edge, higher motion correlation
       Motion_indices = normalize(motion_corr) > 1 & normalize(ROI_props(:,2)) > 0;
       Motion_ICs = ICs(Motion_indices);

       % CSF or internal artery-like features: CSF is hard to catch, so:
       % high freq, clustering, and outbrain
       Phys_indices = sum(normalize(ROI_props(:,3:4) ./ ROI_props(:,1))>1,2) > 0;
       Phys_ICs = ICs(Phys_indices);

       % Susceptibility-like features: high Suscept, high low freq, low BOLD and high
       % freq
       Suscept_indices = normalize(ROI_props(:,5))>1;
       Suscept_ICs = ICs(Suscept_indices);
       
       % MRIart: high frequency and low clustering prop
       MRIArt_indices = normalize(clustering_prop) < -1 & normalize(freq_props(:,3)) > 1;
       MRIArt_ICs = ICs(MRIArt_indices);
       
       % Anything labeled as noise and potential GM is most likely noise
       GM_indices = GM_pot_indices == 1 & Motion_indices == 0 & Phys_indices == 0 ...
           & Suscept_indices == 0 & MRIArt_indices == 0;
       GM_ICs = ICs(GM_indices);

       % Unclassified Noise: Anything not labeled, could be MRI artifact,
       % or a mix, or truly unknown. Might not capture anything
       Unclass_indices = GM_pot_indices == 0 & Motion_indices == 0 & Phys_indices == 0 ...
           & Suscept_indices == 0 & MRIArt_indices == 0;
       Unclass_ICs = ICs(Unclass_indices);

       % Create an easy to read table of classification labels
       classification_names = {'IC Number', 'GM-like', 'Motion/Edge-like', 'Phys-like', 'Suscept-like', 'MRIArt-like', 'Unclassified'};
       classification_array = [ICs', GM_indices, Motion_indices, Phys_indices, Suscept_indices, MRIArt_indices, Unclass_indices];
       classification_table = array2table(classification_array, 'VariableNames', classification_names);

       % Calculate ones that may be close to being GM or not.
       closeGM_indices = (GM_indices == 1 & ((ROI_props(:,1)./max(ROI_props, [], 2)) < 1)) ...
       | (GM_indices == 0 & sum(X_neg < 0.25, 2) == width(X_neg) & clustering_prop > 0.45 & normalize(motion_corr) < 1);
       closeGM_ICs = ICs(closeGM_indices);

       % Make an array and table showing current identification of close ICs
       close_IC_checker = array2table([closeGM_ICs', GM_indices(closeGM_ICs)], 'VariableNames', {'Close IC', 'Labeled as GM Signal?'});
       close_IC_classification = classification_table(closeGM_ICs, :);

       % Overall evaluation averages for each type
       means_array = [mean(X(GM_indices,:),1); mean(X(Motion_indices,:),1); mean(X(Phys_indices,:),1); ...
           mean(X(Suscept_indices,:),1); mean(X(MRIArt_indices,:),1); mean(X(Unclass_indices,:),1)]';
       std_array = [std(X(GM_indices,:),0,1); std(X(Motion_indices,:),0,1); std(X(Phys_indices,:),0,1); ...
           std(X(Suscept_indices,:),0,1); std(X(MRIArt_indices,:),0,1); std(X(Unclass_indices,:),0,1)]';
       averages_table = array2table(means_array, 'VariableNames', classification_names(:,2:end), ...
           'RowNames', factor_names);
       std_table = array2table(std_array, 'VariableNames', classification_names(:,2:end), ...
           'RowNames', factor_names);

       % (3b) Make it easy to compare groups with the normalized table
       eval_names = {'IC Number', 'GM', 'Edge', 'Subependymal', 'CSF', 'Suscept', 'BOLDfreq', 'HighFreq', 'Motion_Corr', 'Clustering_Prop'};
       eval_table = table(ICs', X(:,1), X(:,2), X(:,3), X(:,4), X(:,5), X(:,6), X(:,7), X(:,8), X(:,9), 'VariableNames', eval_names);
       eval_table_close_ICs = eval_table(closeGM_ICs, :);

       % None of the other regions outside of grey matter are leaned to more
       % heavily, and at least 50% of the data significantly clusters
       GM_indices = sum(X_neg < 0, 2) == width(X_neg) & normalize(clustering_prop) > 0.5;
       GM_ICs = ICs(GM_indices);

       % consider: ICs(normalize(ROI_props(:,1)) > 0 & normalize(clustering_prop) > 0 & normalize(Outbrain_prop) < 0 & normalize(motion_corr) < 0)

       % (3c) Create Signal Labels Officially
       final_indices = GM_indices; % this is the variable to change if you want to override some labels
       signal_indices = logical(final_indices);
       signal_ICs = ICs(signal_indices);

       noise_indices = logical(~final_indices);
       noise_ICs = ICs(noise_indices);

       % (3d) Compare before and after selection
       before = sum(feature_array .* IC_exp_var);
       IC_exp_var_signal = IC_exp_var(signal_ICs);
       IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
       after = sum(feature_array(signal_ICs,:) .* IC_exp_var_signal);
       compare_cleaning = array2table([before; after], 'VariableNames', eval_names(2:end), 'RowNames', {'Before', 'After'});

       % (5c) Export variables later use in next steps (e.g., fsl_regfilt)
       writematrix(noise_ICs, 'Noise_dist_ICs.csv')
       writematrix(signal_ICs, 'Signal_dist_ICs.csv')
        
       % and then make noise components if you want to do aggressive denoising like in CONN (regression)
       mixing_matrix = dlmread('./melodic/melodic_mix');
       noise_dist_covariates = mixing_matrix(:, noise_ICs);
       save('CADICA_Noise_dist.mat', 'noise_dist_covariates')

       % save a matrix of relevant variables if you want to examine later!
       save('DecisionVariables.mat', 'ICs', 'classification_table',...
           'X', 'eval_table', 'signal_ICs', 'noise_ICs')
    end
end


% Now, go ahead and compare results and labeling to the FSL outputs in
% Melodic (check report) - you can likely just look at the close ICs. 
% You can open the related files too in fsleyes to
% help see the decisions. We suggest overlaying fullvolICA_adj on the
% mnitemplate


