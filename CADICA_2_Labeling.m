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
cadicafol = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA';
% similar syntax as in the 1_CADICA_MasksandICAs
% subjects = {'102' '103' '105' '106' '107' '108' '109' '112' '114' '115' '116' '117' '121' '122' '125' '126' '128' '129' '130' '132' '134' '135' '136' '138' '139' '140' '142' '143' '144' '145' '146' '147' '149' '153' '156' '157' '158' '160' '161' '164' '165' '168' '169' '171' '172' '173' '174' '176' '178' '179' '181' '184' '185' '186' '187'};
subjects = {'108'};
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

       ROI_props = [GM_prop, Edge_prop, WMCSF_prop, CSF_prop, Suscept_prop]; % we don't need Outbrain or WM
       ROI_props_more = [GM_prop, Outbrain_prop, Edge_prop, WMCSF_prop, CSF_prop, Suscept_prop];
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

       freq_props = [lowfreqIC_prop, BOLDfreqIC_prop, highfreqIC_prop]; %don't need lowfreq
       freq_props_table = array2table(freq_props, 'VariableNames', freqs);

       % (2c) Select freq props based on majority like with GM
       BOLDfreq_dom_indices = freq_props(:,1) == max(freq_props, [], 2);

       % (3) Correlation to movement spikes, using framewise displacement
       confound_place = [confoundext, 'confounds_timeseries'];
       allconfounds = readtable(confound_place);
       confounds_rmsd = allconfounds(:,{'rmsd'});
       rmsd_corr = corr(table2array(confounds_rmsd(2:end,:)), normalize(abs(diff(detrend(ts))), "range"))';

       % (4) Putting it all together!

       % Make a normalized array so it is easier to follow in the future
       feature_array = [ROI_props, freq_props, rmsd_corr.^2, clustering_prop];
       factor_names = {'GM', 'Edge', 'Transmed', 'CSF', 'Suscept', 'low freq', 'BOLD freq', 'high freq', 'Motion corr', 'Clustering prop'};
       feature_table = array2table(feature_array, 'VariableNames', factor_names);
       X = normalize(feature_array, 'medianiqr');

       % GM: GM overlap is largest ROI overlap and higher clustering, BOLD
       % > high freq, and GM higher than edge + CSF
       GM_indices = ROI_props(:,1) == max(ROI_props, [], 2) & clustering_prop > 0.5 ...
           & freq_props(:,2) > freq_props(:,3) & ROI_props(:,1) > (ROI_props(:,2) + ROI_props(:,4));
       GM_ICs = ICs(GM_indices);

       % CSF: high CSF, high high freq, low BOLD and low freq
       CSF_indices = CSF_prop > median(CSF_prop) & highfreqIC_prop > median(highfreqIC_prop)... 
           & lowfreqIC_prop < median(lowfreqIC_prop) & BOLDfreqIC_prop < median(BOLDfreqIC_prop);
       CSF_ICs = ICs(CSF_indices);

       % Susceptibility: high Suscept, high low freq, low BOLD and high
       % freq
       Suscept_indices = Suscept_prop > median(Suscept_prop) & lowfreqIC_prop > median(lowfreqIC_prop)...
           & BOLDfreqIC_prop < median(BOLDfreqIC_prop) & highfreqIC_prop < median(highfreqIC_prop);
       Suscept_ICs = ICs(Suscept_indices);

       % Transmedullary/Subependymal (high in transmedullary, high in BOLD), low in GM,
       % high in clustering. Only one we will say is low in GM, because
       % otherwise very similar.
       Subepe_indices = WMCSF_prop > median(WMCSF_prop) & Edge_prop < median(Edge_prop) ...
            & GM_prop < median(GM_prop) & BOLDfreqIC_prop > median(BOLDfreqIC_prop) & highfreqIC_prop < median(highfreqIC_prop) ...
            & lowfreqIC_prop < median(lowfreqIC_prop) & clustering_prop > median(clustering_prop);
       Subepe_ICs = ICs(Subepe_indices);

       % External Arteries: high edge and BOLD, low highfreq, not CSF or susceptibility
       ExtArt_indices = Edge_prop > median(Edge_prop) & BOLDfreqIC_prop > median(BOLDfreqIC_prop) ...
           & highfreqIC_prop < median(highfreqIC_prop) & CSF_indices == 0 & Suscept_indices == 0;
       ExtArt_ICs = ICs(ExtArt_indices);
       
       % Sinus: not any of the above, but edge and BOLD freq, lower motion
       % corr
        Sinus_indices = Edge_prop > median(Edge_prop) & BOLDfreqIC_prop > median(BOLDfreqIC_prop) ...
           & rmsd_corr < median(rmsd_corr) & ExtArt_indices == 0 & CSF_indices == 0 & Suscept_indices == 0;
        Sinus_ICs = ICs(Sinus_indices);
       
       % Motion: Edge, but not any of the others above
       Motion_indices = Edge_prop > median(Edge_prop) & CSF_prop < median(CSF_prop)...
           & CSF_indices == 0 & Suscept_indices == 0 & ExtArt_indices == 0 & Sinus_indices == 0;
       Motion_ICs = ICs(Motion_indices);
       
       % internal arteries: CSF high, high freq high, but not any of the others
       IntArt_indices = CSF_prop > median(CSF_prop) & CSF_indices == 0 ...
           & Suscept_indices == 0 & highfreqIC_prop > median(highfreqIC_prop) ...
           & Subepe_indices == 0 & ExtArt_indices == 0 & Sinus_indices == 0 & Motion_indices == 0;
       IntArt_ICs = ICs(IntArt_indices);

       % MRIart: high frequency, not any of the others - one may not always
       % find this, so this may be like hidden sinus or something else
       MRIArt_indices = highfreqIC_prop > median(highfreqIC_prop) & BOLDfreqIC_prop < median(BOLDfreqIC_prop) ...
           & lowfreqIC_prop < median(lowfreqIC_prop) & CSF_indices == 0 & IntArt_indices == 0 ...
           & Motion_indices == 0;
       MRIArt_ICs = ICs(MRIArt_indices);

       % Unclassified Noise: Anything not labeled, could be MRI artifact,
       % or a mix, or truly unknown. Might not capture anything
       Unclass_indices = GM_indices == 0 & CSF_indices == 0 & Suscept_indices == 0 ...
           & Subepe_indices == 0 & ExtArt_indices == 0 & Sinus_indices == 0 ...
           & Motion_indices == 0 & IntArt_indices == 0 & MRIArt_indices == 0;
       Unclass_ICs = ICs(Unclass_indices);

       % Create an easy to read table of classification labels
       classification_names = {'IC Number', 'GM-like', 'Motion-like', 'Sinus-like', 'External Artery-like', ...
           'Internal Artery-like', 'CSF-like', 'Subepe-like', 'Suscept-like', 'MRIArt-like', 'Unclassified'};
       classification_array = [ICs', GM_indices, Motion_indices, Sinus_indices, ExtArt_indices, IntArt_indices, ...
           CSF_indices, Subepe_indices, Suscept_indices, MRIArt_indices, Unclass_indices];
       classification_table = array2table(classification_array, 'VariableNames', classification_names);

       % Calculate ones that may be close to being GM or not. We just move
       % the sliding bars a lot, and can examine current GM labels if their
       % GM_prop barely makes the cut. Also if something is labeled both as
       % GM and/or noise.
       closeGM_indices = ((GM_prop + 0.1) > (ROI_props(:,2)+ROI_props(:,4)) & ((GM_prop - 0.1) < (ROI_props(:,2)+ROI_props(:,4)) | GM_indices == 0) ...
           & clustering_prop > 0.4 & (BOLDfreqIC_prop+0.1) > highfreqIC_prop) | (classification_array(:,2) == 1 & sum(classification_array(:,3:end),2) > 0);
       closeGM_ICs = ICs(closeGM_indices);

       % Make an array and table showing current identification of close ICs
       close_IC_checker = array2table([closeGM_ICs', GM_indices(closeGM_ICs)], 'VariableNames', {'Close IC', 'Labeled as GM Signal?'});
       close_IC_classification = classification_table(closeGM_ICs, :);

       % Overall evaluation averages for each type
       means_array = [mean(X(GM_indices,:),1); mean(X(Motion_indices,:),1); mean(X(Sinus_indices,:),1); ...
           mean(X(ExtArt_indices,:),1); mean(X(IntArt_indices,:),1); mean(X(CSF_indices,:),1); mean(X(Subepe_indices,:),1); ...
           mean(X(Suscept_indices,:),1); mean(X(MRIArt_indices,:),1); mean(X(Unclass_indices,:),1)]';
       std_array = [std(X(GM_indices,:),0,1); std(X(Motion_indices,:),0,1); std(X(Sinus_indices,:),0,1); ...
           std(X(ExtArt_indices,:),0,1); std(X(IntArt_indices,:),0,1); std(X(CSF_indices,:),0,1); std(X(Subepe_indices,:),0,1); ...
           std(X(Suscept_indices,:),0,1); std(X(MRIArt_indices,:),0,1); std(X(Unclass_indices,:),0,1)]';
       averages_table = array2table(means_array, 'VariableNames', classification_names(:,2:end), ...
           'RowNames', factor_names);
       std_table = array2table(std_array, 'VariableNames', classification_names(:,2:end), ...
           'RowNames', factor_names);

       % (3b) Make it easy to compare groups with the normalized table
       eval_names = {'IC Number', 'GM', 'Edge', 'Subependymal', 'CSF', 'Suscept', 'low freq' 'BOLDfreq', 'HighFreq', 'Motion_Corr', 'Clustering_Prop'};
       eval_table = table(ICs', X(:,1), X(:,2), X(:,3), X(:,4), X(:,5), X(:,6), X(:,7), X(:,8), X(:,9), X(:,10), 'VariableNames', eval_names);
       eval_table_close_ICs = eval_table(closeGM_ICs, :);

       % (3c) Create Signal Labels Officially
       signal_indices = logical(GM_indices);
       signal_ICs = ICs(signal_indices);

       noise_indices = logical(~GM_indices);
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


