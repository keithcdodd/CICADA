%% Script to be run after 1_CADICA_MasksandICAs_MNI.sh
% This script increases efficiency for selecting signal vs noise ICs. It
% will make strong educated guesses on noise vs signal labels, but it is ideally used to help inform
% the user what IC's should be double-checked

%% How To Use:
% 1. Make sure starting variables are all appropriate. Then run it. You can
% run it for many subjects/sessions at a time if you would like.
% 2. clearvars and then load the DecisionVariables.mat for the
% subject/session you want to update (this step is not necessary if you
% only ran this script on one subject/session/run).
% 2. Open up the variables IC_checker_updated table, and
% potential_IC_classification variables in Matlab.
% 3. Open up the corresonding MELODIC report (html), and the anatfile and
% melodic_IC.nii.gz overlayed in FSLEYES
% 4. Update signal label of each IC in the IC_checker_updated_table, use
% the information in the IC_checker_updated_table, the MELODIC report, and
% the overlay in FSLEYES to help inform your decision. 
% 5. After updating the signal labels in IC_checker_updated_table, evaluate
% the indicated section of this script in the command window again (see 2nd
% half of this script).
% 6. You are done. Proceed to the next script.

%% General Method of This Script:
% This script uses the proportion of GM overlap, the ratio of BOLD power to
% high frequency power, and the smoothness of each IC to determine a rank
% order probability of true signal vs noise labeling ICs. IC's that should be
% checked then are given by the top ranked ICs. The script uses
% measurements of known potential noise variables to decide when to stop
% adding ICs to the IC checker. Therefore, this script not only auto-labels
% potential signal vs noise ICs, but it also greatly decreases the time to
% evaluate ICs manually by significantly decreasing the number of ICs needed to
% check, and by providing more helpful information to evaluate the ICs.

clearvars

%%%%% Everything Else
% cadicafol = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA_Updated';
cadicafol='/Volumes/VectoTec_VectoTech_Media_Rapid/DMXBA/CADICA_Updated';
% similar syntax as in the 1_CADICA_MasksandICAs
% subjects = {'102' '103' '105' '106' '107' '109' '112' '114' '115' '116' '117' '121' '122' '125' '126' '128' '129' '130' '132' '134' '135' '136' '138' '139' '140' '142' '143' '144' '145' '146' '147' '149' '153' '156' '157' '158' '160' '161' '164' '165' '168' '169' '171' '172' '173' '174' '176' '178' '179' '181' '184' '185' '186' '187'};
subjects = {'067'};
% subjects = {'102', '103'};
sessions = {'03'};
tasks = {'rest'}; % or 'rest', this needs to match associated folder
task_hrfs = struct([]); % load([cadicafol, '/designmatrix2.mat']); % make struct([]) if resting state
task_names = {}; % {'Highcal', 'Lowcal', 'Nonfood'}; % if there are task_hrfs, list these in the same order, otherwise if resting, put {};
numvolumes = 300; % 172; % 300 % how many TRs/samples, will be different depending on task
TR = 2; % in seconds
ROIs = {'GM' 'Edge' 'Transmedullary' 'CSF', 'Susceptibility', 'Outbrain'};
freqs = {'lowfreq', 'BOLDfreq', 'higherfreq'};
% Need to know whether or not to take frequency power into account (was it
% already bandpassed?)
bandpassed = 0; % 1 if it was bandpassed before MELODIC, 0 if it was not
funcfilename = 'funcfile';
%funcfilename = 's_bp_funcfile';
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
        for l = 1:length(tasks)
            cd(cadicafol)
            cd(currsubjfol)
            cd(currsessfol)
            currtask = tasks{l};
            cd(currtask)

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

            % smoothness (min cluster size) -> might be useful. Make
            % inverse so that a higher number does mean more smoothness
            smoothness_vals = 1 ./ dlmread('clustering/smoothness.txt');
            minclustersize = dlmread('clustering/clustersizes_pos.txt');

            % Outbrain with CSF
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

            % White Matter
            WM_ICmean = dlmread('ROIcalcs/WM_ICmean.txt');
            WM_ICnumvoxels = dlmread('ROIcalcs/WM_ICnumvoxels.txt');
            WM_ICsum = WM_ICmean .* WM_ICnumvoxels;

            % White Matter CSF Boundary (Subependymal)
            WMCSF_ICmean = dlmread('ROIcalcs/Subepe_ICmean.txt');
            WMCSF_ICnumvoxels = dlmread('ROIcalcs/Subepe_ICnumvoxels.txt');
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
            LowfreqIC = zeros(length(ICs), 1);
            HighfreqIC = zeros(length(ICs), 1);

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
                LowfreqIC(i)=f_lowfreq;
                HighfreqIC(i)=f_highfreq;
            end

            % (2b) Calculate actual measured proportions.
            LowfreqIC_prop = LowfreqIC ./ (LowfreqIC+BOLDfreqIC+HighfreqIC);
            BOLDfreqIC_prop = BOLDfreqIC ./ (LowfreqIC+BOLDfreqIC+HighfreqIC);
            HighfreqIC_prop = HighfreqIC ./ (LowfreqIC+BOLDfreqIC+HighfreqIC);

            freq_props = [LowfreqIC_prop, BOLDfreqIC_prop, HighfreqIC_prop];
            freq_props_table = array2table(freq_props, 'VariableNames', freqs);

            % (3) Correlations: to dvars, FD, tasks (if not resting state), and relevant noise areas:
            confound_place = 'confounds_timeseries.csv';
            allconfounds = readtable(confound_place);
            confounds_9param = table2array(allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'white_matter', 'csf', 'global_signal'}));
            confounds_8param = table2array(allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'white_matter', 'csf'}));
            confounds_dvars = table2array(allconfounds(:,{'dvars'}));
            confounds_fd = table2array(allconfounds(:,{'framewise_displacement'}));

            fd_corr = corr(confounds_fd(2:end), abs(diff(detrend(ts)))).^2';
            dvars_corr = corr(confounds_dvars(2:end), abs(diff(detrend(ts)))).^2';

            % if there are task hrf functions (which could have been
            % created in spm and then saved as a .mat in CADICA folder
            % ahead of time) then also look at correlations to those, since
            % those can help indicate ICs that are potentially relevant
            task_corrs = [];
            if ~isempty(fieldnames(task_hrfs))
                task_hrfs_table = struct2table(task_hrfs);
                task_hrfs_table.Combined = sum(task_hrfs_table{:,:},2);
                task_names = [task_names, 'Combined'];
                task_corrs = zeros(length(ICs), length(fieldnames(task_hrfs)));
                for m = 1:width(task_hrfs_table)
                    task_corrs(:,m) = corr(table2array(task_hrfs_table(:,m)), ts).^2;
                end
            end

            % save standard 9Param and 8Param for comparison or use later:
            writematrix(confounds_9param, ['9p_regressors.txt'], 'Delimiter', ' ')
            writematrix(confounds_8param, ['8p_regressors.txt'], 'Delimiter', ' ')

            % (4) Putting it all together!

            % Make a normalized array so it is easier to follow in the future
            feature_array = [ROI_props, freq_props, dvars_corr, fd_corr, clustering_prop, task_corrs];
            factor_names = ['GM', 'Edge', 'Transmed', 'CSF', 'Suscept', 'Outbrain', 'low freq', 'BOLD freq', 'high freq', 'Spike corr', 'Move corr', 'Clustering prop', task_names];
            feature_table = array2table(feature_array, 'VariableNames', factor_names);
            feature_array_norm = normalize(feature_array);
            feature_table_norm = array2table(feature_array_norm, 'VariableNames', factor_names);


            % Types of Noise to remove - use a for loop to find ones that are
            % close to being classified as signal or noise
            noise_labelling = zeros(length(ICs), 3);

            % let's try with kmeans clustering for each, incorporate signal
            % measures from before too

            % Edge
            curr_prop = Edge_prop ./ (Edge_prop + GM_prop);
            [prop_idx, prop_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            HighEdge_indices = prop_idx == find(prop_C == max(prop_C));
            LowEdge_indices = prop_idx == find(prop_C == min(prop_C));
            Edge_ICs = ICs(HighEdge_indices);

            % Subependymal (WMCSF)
            curr_prop = WMCSF_prop ./ (GM_prop+WMCSF_prop);
            [prop_idx, prop_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            HighWMCSF_indices = prop_idx == find(prop_C == max(prop_C));
            LowWMCSF_indices = prop_idx == find(prop_C == min(prop_C));
            WMCSF_ICs = ICs(HighWMCSF_indices);

            % CSF
            curr_prop = CSF_prop ./ (GM_prop+CSF_prop);
            [prop_idx, prop_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            HighCSF_indices = prop_idx == find(prop_C == max(prop_C));
            LowCSF_indices = prop_idx == find(prop_C == min(prop_C));
            CSF_ICs = ICs(HighCSF_indices);

            % Outbrain
            curr_prop = Outbrain_prop;
            [prop_idx, prop_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            HighOutbrain_indices = prop_idx == find(prop_C == max(prop_C));
            LowOutbrain_indices = prop_idx == find(prop_C == min(prop_C));
            Outbrain_ICs = ICs(HighOutbrain_indices);

            % Low Frequency
            curr_prop = LowfreqIC_prop ./ BOLDfreqIC_prop;
            [Lowfreq_idx, Lowfreq_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            HighLowfreq_indices = Lowfreq_idx == find(Lowfreq_C == max(Lowfreq_C));
            LowLowfreq_indices = Lowfreq_idx == find(Lowfreq_C == min(Lowfreq_C));
            Lowfreq_ICs = ICs(HighLowfreq_indices);

            % High Frequency
            curr_prop = HighfreqIC_prop ./ BOLDfreqIC_prop;
            [Highfreq_idx, Highfreq_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            HighHighfreq_indices = Highfreq_idx == find(Highfreq_C == max(Highfreq_C));
            LowHighfreq_indices = Highfreq_idx == find(Highfreq_C == min(Highfreq_C));
            Highfreq_ICs = ICs(HighHighfreq_indices);

            % dvars
            curr_prop = dvars_corr;
            [dvarscorr_idx, dvarscorr_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            Highdvarscorr_indices = dvarscorr_idx == find(dvarscorr_C == max(dvarscorr_C));
            Lowdvarscorr_indices = dvarscorr_idx == find(dvarscorr_C == min(dvarscorr_C));
            dvarscorr_ICs = ICs(Highdvarscorr_indices);

            % framewise displacement
            curr_prop = fd_corr;
            [fdcorr_idx, fdcorr_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            Highfdcorr_indices = fdcorr_idx == find(fdcorr_C == max(fdcorr_C));
            Lowfdcorr_indices = fdcorr_idx == find(fdcorr_C == min(fdcorr_C));
            fdcorr_ICs = ICs(Highfdcorr_indices);

            % clustering prop
            curr_prop = clustering_prop;
            [Cluster_idx, Cluster_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            LowCluster_indices = Cluster_idx == find(Cluster_C == min(Cluster_C));
            LowCluster_ICs = ICs(LowCluster_indices);
            HighCluster_indices = Cluster_idx == find(Cluster_C == max(Cluster_C));
            HighCluster_ICs = ICs(HighCluster_indices);

            % smoothness: Because we inverted it, higher numbers are more
            % smooth, which is convenent
            curr_prop = smoothness_vals;
            [Smoothness_idx, Smoothness_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            LowSmoothness_indices = Smoothness_idx == find(Smoothness_C == min(Smoothness_C));
            LowSmoothness_ICs = ICs(LowSmoothness_indices);
            HighSmoothness_indices = Smoothness_idx == find(Smoothness_C == max(Smoothness_C));
            HighSmoothness_ICs = ICs(HighSmoothness_indices);

            % GM
            curr_prop = GM_prop;
            [prop_idx, prop_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            LowGM_indices = prop_idx == find(prop_C == min(prop_C));
            HighGM_indices = prop_idx == find(prop_C == max(prop_C));
            LowGM_ICs = ICs(LowGM_indices);
            HighGM_ICs = ICs(HighGM_indices);

            % BOLD frequency
            curr_prop = BOLDfreqIC_prop;
            [LowBOLDfreq_idx, LowBOLDfreq_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
            LowBOLDfreq_indices = LowBOLDfreq_idx == find(LowBOLDfreq_C == min(LowBOLDfreq_C));
            LowBOLDfreq_ICs = ICs(LowBOLDfreq_indices);
            HighBOLDfreq_indices = LowBOLDfreq_idx == find(LowBOLDfreq_C == max(LowBOLDfreq_C));
            HighBOLDfreq_ICs = ICs(HighBOLDfreq_indices);

            % task correlated
            Lowtaskcorr_indices = []; Hightaskcorr_indices = [];
            Hightaskcorr_ICs = []; Lowtaskcorr_ICs = [];
            if ~isempty(fieldnames(task_hrfs))
                Hightaskcorr_indices = zeros(length(ICs), length(fieldnames(task_hrfs)));
                Lowtaskcorr_indices = Hightaskcorr_indices;
                for m = 1:width(task_hrfs_table)
                    curr_prop = task_corrs(:,m);
                    [taskcorr_idx, taskcorr_C] = kmeans(curr_prop, 3, 'Start', [min(curr_prop); median(curr_prop); max(curr_prop)]);
                    Hightaskcorr_indices(:,m) = taskcorr_idx == find(taskcorr_C == max(taskcorr_C));
                    Lowtaskcorr_indices(:,m) = taskcorr_idx == find(taskcorr_C == min(taskcorr_C));
                    Hightaskcorr_ICs = [Hightaskcorr_ICs, ICs(logical(Hightaskcorr_indices(:,m)))];
                    Lowtaskcorr_ICs = [Lowtaskcorr_ICs, ICs(logical(Lowtaskcorr_indices(:,m)))];
                end

                Hightaskcorr_ICs = unique(Hightaskcorr_ICs);
                Lowtaskcorr_ICs = unique(Lowtaskcorr_ICs);
            end

            % Now calculate ICs that are almost certainly High Noise or High
            % Signal
            if bandpassed == 1
                % if it is bandpassed, low and high frequency power are not
                % reliable indicators
                % high signal
                HighSig_indices = logical(HighGM_indices .* HighCluster_indices .* ...
                    LowEdge_indices .* LowWMCSF_indices .* LowCSF_indices .* ...
                    LowOutbrain_indices);
                HighSig_ICs = ICs(HighSig_indices);

                % high noise
                Highnoise_indices = logical(HighEdge_indices + HighWMCSF_indices + HighCSF_indices + HighOutbrain_indices + ...
                    Highdvarscorr_indices + Highfdcorr_indices + LowCluster_indices);
                Highnoise_ICs = ICs(Highnoise_indices);
            else
                % if it is not bandpassed, then take low and high frequency
                % content into account for signal and noise
                % high signal
                HighSig_indices = logical(HighGM_indices .* HighCluster_indices .* HighBOLDfreq_indices .* ...
                    LowEdge_indices .* LowWMCSF_indices .* LowCSF_indices .* ...
                    LowOutbrain_indices .* LowLowfreq_indices .* LowHighfreq_indices);
                HighSig_ICs = ICs(HighSig_indices);

                % high noise
                Highnoise_indices = logical(HighEdge_indices + HighWMCSF_indices + HighCSF_indices + HighOutbrain_indices + LowGM_indices...
                    + Highdvarscorr_indices + Highfdcorr_indices + HighHighfreq_indices + HighLowfreq_indices + LowBOLDfreq_indices + LowCluster_indices);
                Highnoise_ICs = ICs(Highnoise_indices);
            end

            % OK, now we can start to narrow down our search! Good ICs are
            % determined by GM overlap (compared to other regions),
            % Boldfrequency power (compared to high frequency), and
            % smoothness. Weight GM overlap heavier than smoothness or
            % frequency (because noise could be in BOLD and smooth, like
            % WMCSF signal)
            [signal_val, signal_idx] = sort(normalize(smoothness_vals, "range") .* normalize(GM_prop, "range").^2 .* normalize(BOLDfreqIC_prop./HighfreqIC_prop, "range"), 'descend');

            % OK, step through, and if 3 in a row are high noise, then stop
            % looking:
            tolerance = 3;
            g = 1;
            noisecount = 0;
            while noisecount < tolerance
                if ismember(signal_idx(g), Highnoise_ICs) == 1
                    noisecount = noisecount + 1;
                else
                    noisecount = 0;
                end
                g = g + 1; % g will be how far along the signal ICs you go
            end
            g = g - 1;

            potential_signal_ICs = signal_idx(1:g);
            potential_signal_indices = zeros(length(ICs),1);
            potential_signal_indices(potential_signal_ICs) = 1;
            potential_signal_indices = logical(potential_signal_indices);

            % if they are labeled as highnoise, switch them to that label
            signal_indices = logical(potential_signal_indices .* ~Highnoise_indices);
            noise_indices = ~signal_indices;
            signal_ICs = ICs(signal_indices);
            Noise_ICs = ICs(noise_indices);

            signal_checker_table = array2table(potential_signal_ICs, 'VariableNames', {'Potential ICs'});
            signal_checker_table.SignalLabel = signal_indices(potential_signal_ICs);
            signal_checker_table.HighSignalLabel = HighSig_indices(potential_signal_ICs);
            signal_checker_table.HighNoiseLabel = Highnoise_indices(potential_signal_ICs);
            if ~isempty(fieldnames(task_hrfs))
                signal_checker_table.TaskCorrelated = double(sum(Hightaskcorr_indices(potential_signal_ICs,:),2) > 0);
            end        

            % Create an easy to read table of classification labels
            factors = [signal_indices, HighEdge_indices, HighWMCSF_indices, HighCSF_indices, HighOutbrain_indices, ...
                HighLowfreq_indices, HighHighfreq_indices, Highdvarscorr_indices, Highfdcorr_indices, LowCluster_indices, Hightaskcorr_indices];
            classification_names = ['IC Number', 'GM-like', 'Edge-like', 'WMCSF-like', 'CSF-like', ...
                'Outbrain-like', 'Low Freq', 'High Freq', 'Spike Correlated', 'Movement Correlated', 'Not Clustered', task_names];
            classification_array = [ICs', factors];
            classification_table = array2table(classification_array, 'VariableNames', classification_names);

            potential_IC_classification = classification_table(sort(potential_signal_ICs)', :);


            % Overall evaluation averages for each type (before final editing)
            for n = 1:size(factors,2)
                means_array(:,n) = mean(feature_array(logical(factors(:,n)),:),1);
                stds_array(:,n) = std(feature_array(logical(factors(:,n)),:),1);
            end
            averages_table = array2table(means_array, 'VariableNames', classification_names(:,2:end), ...
                'RowNames', factor_names);
            std_table = array2table(stds_array, 'VariableNames', classification_names(:,2:end), ...
                'RowNames', factor_names);


            IC_checker_updated_table = signal_checker_table;


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % After updating the signal column labels in
            % IC_checker_updated_table, re-evaluate all code below in the
            % command window. Alternatively, if you ran this script on
            % multiple subjects/sessions at once, you will want to go into
            % each subject folder, load the DecisionVariables.mat,
            % update the IC_checker_updated table, and THEN re-evaluate all
            % code below in the command window (not including the ends of
            % the for loops)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            %%%%% Update IC_checker_updated_table variable above before
            %%%%% re-evaluating everything below
            close_IC_checker_updated_array = table2array(IC_checker_updated_table);
            noise_indices(close_IC_checker_updated_array(close_IC_checker_updated_array(:,2) == 0,1)) = 1;
            noise_indices(close_IC_checker_updated_array(close_IC_checker_updated_array(:,2) == 1,1)) = 0;
            noise_indices = logical(noise_indices);

            % Create final noise and signal ICs and indices, if GM_prop>0.5,
            % label as signal regardless
            % check close ICs and update noise indices
            noise_ICs = ICs(noise_indices);
            signal_indices = logical(~noise_indices);
            signal_ICs = ICs(signal_indices);


            % (3b) Make it easy to compare groups
            close_ICs = int32(close_IC_checker_updated_array(:,1));
            eval_names = ['IC Number', 'GM', 'Edge', 'Subependymal', 'CSF', 'Suscept', 'Outbrain', 'low freq' 'BOLDfreq', 'HighFreq', 'Spike_Corr', 'Move _Corr', 'Clustering_Prop', task_names];
            eval_table = array2table([ICs', feature_array], 'VariableNames', eval_names);
            eval_table_close_ICs = eval_table(close_ICs, :);

            % (3d) Compare before and after selection
            before = sum(feature_array .* IC_exp_var);
            IC_exp_var_signal = IC_exp_var(signal_ICs);
            IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
            after = sum(feature_array(signal_ICs,:) .* IC_exp_var_signal);
            compare_cleaning = array2table([before', after'], 'RowNames', eval_names(2:end), 'VariableNames', {'Before', 'After'});

            % estimate effective degrees of freedom (somewhat) from CADICA alone via proportions:
            percent_variance_kept = sum(IC_exp_var(signal_ICs)); % probably want > 0.1 (10%)

            % Create images to review signal and noise components from
            % probabilities ICs
            noise_prob = niftiread('./melodic/ICprobabilities.nii.gz');
            noise_prob_info = niftiinfo('./melodic/ICprobabilities.nii.gz');
            signal_prob = noise_prob; signal_prob_info = noise_prob_info;

            % grab noise probabilities, update header, and write, same for
            % signal
            noise_prob = noise_prob(:,:,:,noise_ICs);
            noise_prob_info.ImageSize = size(noise_prob);
            noise_prob_info.DisplayIntensityRange = [0.9 1];
            signal_prob = signal_prob(:,:,:, signal_ICs);
            signal_prob_info = noise_prob_info;
            signal_prob_info.ImageSize = size(signal_prob);

            niftiwrite(noise_prob, 'NoiseICs', noise_prob_info, 'Compressed', true)
            niftiwrite(signal_prob, 'SignalICs', signal_prob_info, 'Compressed', true)

            % Now binarize and compress noise and signal
            noise_prob(noise_prob>0.899) = 1; noise_prob(noise_prob<0.9) = 0;
            noise_prob = sum(noise_prob, 4); noise_prob(noise_prob>0.99) = 1;

            % grab, binarize and sum signal_probs
            signal_prob(signal_prob>0.899) = 1; signal_prob(signal_prob<0.9) = 0;
            signal_prob = sum(signal_prob, 4); signal_prob(signal_prob>0.99) = 1;

            % Update headers for the 3D images
            noise_prob_info.ImageSize = noise_prob_info.ImageSize(1:end-1);
            noise_prob_info.PixelDimensions = noise_prob_info.PixelDimensions(1:end-1);
            signal_prob_info = noise_prob_info;

            % write out the masks
            niftiwrite(noise_prob, 'NoiseICOverlap', noise_prob_info, 'Compressed', true)
            niftiwrite(signal_prob, 'SignalICOverlap', signal_prob_info, 'Compressed', true)

            % (5c) Export variables later use in next steps (e.g., fsl_regfilt)
            % Auto prefix tells you that it is using the automatic best
            % guess, and not with manual checking. This should be good, but
            % not perfect.
            writematrix(noise_ICs, 'Auto_Noise_dist_ICs.csv')
            writematrix(signal_ICs, 'Auto_Signal_dist_ICs.csv')
            writematrix(noise_indices, 'Auto_Noise_dist_indices.csv')
            writematrix(signal_indices, 'Auto_Signal_dist_indices.csv')

            % and then make noise components if you want to do aggressive denoising like in CONN (regression)
            mixing_matrix = dlmread('./melodic/melodic_mix');
            Auto_noise_dist_covariates = mixing_matrix(:, noise_ICs);
            save('Auto_CADICA_Noise_dist.mat', 'Auto_noise_dist_covariates')

            % save a matrix of relevant variables if you want to examine later!
            save('DecisionVariables_Auto.mat', 'ICs', 'classification_table', 'potential_IC_classification', ...
                'eval_table', 'signal_ICs', 'noise_ICs', 'signal_checker_table', 'IC_checker_updated_table', ...
                'feature_array', 'feature_table', 'IC_exp_var', 'compare_cleaning', 'close_ICs', 'noise_indices', 'HighSig_ICs', ...
                'Highnoise_ICs', 'averages_table', 'std_table', 'classification_table', 'percent_variance_kept', 'Hightaskcorr_ICs', 'feature_table_norm')
        end
    end
end