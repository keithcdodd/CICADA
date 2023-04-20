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
TR = 2; % in seconds
numvolumes = 300; % how many TRs/samples
cadicafol = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA';
% similar syntax as in the 1_CADICA_MasksandICAs
% subjects = {'102' '103' '105' '106' '107' '108' '109' '112' '114' '115' '116' '117' '121' '122' '125' '126' '128' '129' '130' '132' '134' '135' '136' '138' '139' '140' '142' '143' '144' '145' '146' '147' '149' '153' '156' '157' '158' '160' '161' '164' '165' '168' '169' '171' '172' '173' '174' '176' '178' '179' '181' '184' '185' '186' '187'};
subjects = {'108'};
sessions = {'01'};
ROIs = {'GM' 'Edge' 'Transmedullary' 'CSF', 'Susceptibility'};
freqs = {'lowfreq', 'BOLDfreq', 'higherfreq'};
% always include at least GM, WM, CSF, and outbrain
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
       
       % (1) Check Spatial Map Overlap with ROIs
       % (1a) Approximate overlap with mean and numvoxels

       % outbrain
       Outbrain_ICmean = dlmread('./Outbrain_ICmean.txt');
       Outbrain_ICnumvoxels = dlmread('./Outbrain_ICnumvoxels.txt');
       Outbrain_ICsum = Outbrain_ICmean .* Outbrain_ICnumvoxels;

       % Edge
       Edge_ICmean = dlmread('./Edge_ICmean.txt');
       Edge_ICnumvoxels = dlmread('./Edge_ICnumvoxels.txt');
       Edge_ICsum = Edge_ICmean .* Edge_ICnumvoxels;

       % Grey Matter
       GM_ICmean = dlmread('./GM_ICmean.txt');
       GM_ICnumvoxels = dlmread('./GM_ICnumvoxels.txt');
       GM_ICsum = GM_ICmean .* GM_ICnumvoxels;

       % White Matter
       WM_ICmean = dlmread('./WM_ICmean.txt');
       WM_ICnumvoxels = dlmread('./WM_ICnumvoxels.txt');
       WM_ICsum = WM_ICmean .* WM_ICnumvoxels;

       % White Matter CSF Boundary (for subependymal)
       WMCSF_ICmean = dlmread('./WMCSF_ICmean.txt');
       WMCSF_ICnumvoxels = dlmread('./WMCSF_ICnumvoxels.txt');
       WMCSF_ICsum = WMCSF_ICmean .* WMCSF_ICnumvoxels;

       % CSF
       CSF_ICmean = dlmread('./CSF_ICmean.txt');
       CSF_ICnumvoxels = dlmread('./CSF_ICnumvoxels.txt');
       CSF_ICsum = CSF_ICmean .* CSF_ICnumvoxels;

       % Susceptibility 
       Suscept_ICmean = dlmread('./Suscept_ICmean.txt');
       Suscept_ICnumvoxels = dlmread('./Suscept_ICnumvoxels.txt');
       Suscept_ICsum = Suscept_ICmean .* Suscept_ICnumvoxels;

       % How many ICs are we examining?
       ICs = 1:length(GM_ICmean);

       % (1b) Calculate proportion for each brain region for each IC
       Outbrain_prop = Outbrain_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);
       Edge_prop = Edge_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);
       GM_prop = GM_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);
       WM_prop = WM_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);
       WMCSF_prop = WMCSF_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);
       CSF_prop = CSF_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);
       Suscept_prop = Suscept_ICsum ./ (Outbrain_ICsum + GM_ICsum + WM_ICsum + Suscept_ICsum);

       ROI_props = [GM_prop, Edge_prop, WMCSF_prop, CSF_prop, Suscept_prop]; % we don't need Outbrain or WM
       ROI_props(isinf(ROI_props)|isnan(ROI_props)) = 0; % replace nonsensicals
       ROI_props_table = array2table(ROI_props, 'VariableNames', ROIs);

       % (2) Check Power Frequency Analysis in low, BOLD (0.01-0.1), and high
       % frequencies
       % (2a) Grab time series and calculate power and proportions in each
       % frequency range of interest
       N = T/dt;
       F = 1/dt;
       df = 1/T;
       N_freq = N/2 + 1;
       lower_phys_cutoff = round(0.01 / df) + 1; % start at freq 0 at position 1
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
       allconfounds = readtable('confounds_timeseries');
       confounds_dvars = allconfounds(:,{'std_dvars'});
       dvars_corr = corr(table2array(confounds_dvars(2:end,:)), normalize(abs(diff(detrend(ts)))))'; % sign matters for std_dvars, not ts (magnitude of change matters)

       % (4) Putting it all together!
       X = normalize([ROI_props, freq_props, dvars_corr.^2]);


       % GM: favors GM and BOLD. If it super favors GM, it is likely GM too
       GM_indices_1 = (X(:,1) > 0 & X(:,7) > 0);
       GM_indices_2 = (X(:,1) > 1 & X(:,7) > -1);
       GM_indices = logical(GM_indices_1 + GM_indices_2);
       GM_ICs = ICs(GM_indices);

       % Motion: favors edge and motion corr
       motion_indices = (X(:,2) > 0 & X(:,9) > 0);
       motion_ICs = ICs(motion_indices);

       % Deep veins: favors subependymal, low edge, high BOLD, low GM
       subepe_indices = (X(:,3) > 0 & X(:,2) < 0 & X(:,7) > 0 & X(:,1) < 0);
       subepe_ICs = ICs(subepe_indices);

       % CSF: favors CSF and high frequency
       CSF_indices = (X(:,4) > 0 & X(:,8) > 0);
       CSF_ICs = ICs(CSF_indices);

       % External veins & arteries: favors edge, BOLD, low GM
       sinus_indices = (X(:,2) > 0 & X(:,7) > 0 & X(:,1) < 0);
       sinus_ICs = ICs(sinus_indices);

       % Deep Arteries: favors CSF, has some mix of low and high freq
       arteries_indices = (X(:,4) > 0 & abs(X(:,6)) < 1 & abs(X(:,7)) < 1 & abs(X(:,8)) < 1);
       arteries_ICs = ICs(arteries_indices);

       % susceptibility: favors susc and low frequency
       susc_indices = (X(:,5) > 0 & X(:,6) > 0);
       susc_ICs = ICs(susc_indices);
       
       % MRI artifact: entirely high freq, no CSF preference, no motion correlation
       mriart_indices = (X(:,8) > 1 & X(:,4) < 0 & X(:,9) < 0);
       mriart_ICs = ICs(mriart_indices);

       % Unclassified Noise: Calculate from whats left
       noiselabel_indices = logical(motion_indices + CSF_indices + subepe_indices + ...
           sinus_indices + arteries_indices + susc_indices + mriart_indices);
       unclassified_noise_indices = (~noiselabel_indices - GM_indices) == 1;
       unclassified_noise_ICs = ICs(unclassified_noise_indices);
       noise_label_indices = logical(noiselabel_indices + unclassified_noise_indices);

       % Calculate ones that were labeled as GM and as potentially
       % something else
       double_labeled_indices = (GM_indices + noiselabel_indices) > 1;
       close_GM_cutoff_indices = X(:,1) > -0.5 & X(:,1) < 0 & X(:,7) > 0;
       close_BOLD_cutoff_indices = X(:,7) > -0.5 & X(:,7) < 0 & X(:,1) > 0;
       close_indices = logical(double_labeled_indices + close_GM_cutoff_indices + close_BOLD_cutoff_indices);
       close_ICs = ICs(close_indices);

       % Make a table showing current identification of close ICs
       close_IC_checker = array2table([close_ICs', GM_indices(close_ICs)], 'VariableNames', {'Close IC', 'Labeled as Signal?'});


       % (3a) Create an easy to read table of classification labels
       classification_names = {'GM-like', 'Motion-like', 'Sinus-like', 'Artery-like', ...
           'CSF-like', 'Transmed-like', 'Suscept-like', 'MRIArt-like', 'Unclassified'};
       classification_array = [GM_indices, motion_indices, sinus_indices, arteries_indices, ...
           CSF_indices, subepe_indices, susc_indices, mriart_indices, unclassified_noise_indices];
       classification_table = array2table(classification_array, 'VariableNames', classification_names);

       % Overall evaluation averages for each type
       factor_names = {'GM', 'Edge', 'Transmed', 'CSF', 'Suscept', 'low freq', 'BOLD freq', 'high freq', 'Motion corr'};
       means_array = [mean(X(GM_indices,:),1); mean(X(motion_indices,:),1); mean(X(sinus_indices,:),1); ...
           mean(X(arteries_indices,:),1); mean(X(CSF_indices,:),1); mean(X(subepi_indices,:),1); ...
           mean(X(susc_indices,:),1); mean(X(mriart_indices,:),1); mean(X(unclassified_noise_indices,:),1)]';
       averages_table = array2table(means_array, 'VariableNames', classification_names, ...
           'RowNames', factor_names);

       % (3b) Make it easy to compare groups with the normalized table
       eval_names = {'GM', 'Edge', 'Transmed', 'CSF', 'Suscept', 'low freq' 'BOLDfreq', 'HighFreq', 'DVARS_Corr'};
       eval_table = table(X(:,1), X(:,2), X(:,3), X(:,4), X(:,5), X(:,6), X(:,7), X(:,8), X(:,9), 'VariableNames', eval_names);

       % (3c) Create Signal Labels Officially
       signal_indices = logical(GM_indices);
       signal_ICs = ICs(signal_indices);

       noise_indices = logical(~GM_indices);
       noise_ICs = ICs(noise_indices);

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


