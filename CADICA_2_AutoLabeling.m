function CADICA_2_AutoLabeling(cadicafol, subject_id, session_id, taskname, numvolumes, repetitiontime, tolerance, task_events_file)
% Following CADICA_1 script on an fmriprep data set
% example: autoICAlabeling('/path/CADICA_MNI152Lin6Asym_2mm', '102', '01', 'rest_run-01', 300, 2, 4)
% tolerance is how many noise components to run into in a row before it
% stops looking at the next IC (higher tolerance number gives you more ICs
% to examine. Default is 4)
% condition names and condition_hrfs are optional, but helpful if not
% resting state
% task_events_file: Should be the events.tsv (found in bids_data dir) that
% should be included in your bids formatted data. This is not necessary,
% just slightly potentially helpful.
% relevant mat file e.g., load(designmatrix.mat)

% General Method of This Script:s
% This script uses the proportion of GM overlap, the ratio of BOLD power to
% high frequency power, and the smoothness of each IC to determine a rank
% order probability of true signal vs noise labeling ICs. IC's that should be
% checked then are given by the top ranked ICs. The script uses
% measurements of known potential noise variables to decide when to stop
% adding ICs to the IC checker. Therefore, this script not only auto-labels
% potential signal vs noise ICs, but it also greatly decreases the time to
% evaluate ICs manually by significantly decreasing the number of ICs needed to
% check, and by providing more helpful information to evaluate the ICs.

% How to use: call this function with the proper inputs. This you could
% create your own short matlab script and use for loops, for example, to
% run through many subjects/sessions/tasks. Run this function on subjects
% after CADICA_1 has been run, but before trying to run 2B or later
% functions/scripts.


% initialize
Tables = struct;
Data.tr = repetitiontime;
Data.numvolumes = numvolumes;

T = repetitiontime .* numvolumes;
TR = repetitiontime;
dt = repetitiontime;
funcfilename = 'funcfile';
ROIs = {'GM', 'WM', 'Edge', 'Subependymal', 'CSF', 'Susceptibility', 'OutbrainOnly'};
Networks = {'MedialVisual', 'SensoryMotor', 'DorsalAttention', 'VentralAttention', 'FrontoParietal', 'DefaultModeNetwork', 'Subcortical'};
freqs = {'lowfreq', 'BOLDfreq', 'higherfreq'};
hrf_conditions_table = table([]);

% check optional inputs
if ~exist('task_events_file', 'var')
    task_events_file='';
end

if ~exist('tolerance', 'var')
    tolerance = 4;
end

hrf_plot_norm = spm_hrf(TR); % initialize an hrf function for potential correlations
% should centerpad hrf to have same length
hrf_padded = [hrf_plot_norm', zeros(1,numvolumes-length(hrf_plot_norm))]; % to give relevant resolution of full sampling
Data.HRF_general.hrf_padded = hrf_padded;
if ~strcmp(task_events_file, '')
    if isfile(task_events_file)
        task_events_table = readtable(task_events_file, 'FileType', 'delimitedtext');
        % remove 'baseline' if it exists (do not want to model a baseline)
        baseline_exists = find(strcmp(task_events_table{:,3}, 'baseline'));
        if ~isempty(baseline_exists)
            task_events_table(baseline_exists, :) = [];
        end
        
    else
        fprintf(['Cannont find ', task_events_file, '. Running without task condition information.\n'])
        task_events_file = '';
    return 
    end
end


% Now do things!
cd(cadicafol)
currsubjfol = ['sub-', subject_id];
cd(currsubjfol)
currsessfol = ['ses-', session_id];
cd(currsessfol)
cd(taskname)

% remove old ic_auto_selection folder if it exists
% if isfolder('./ic_auto_selection')
%     rmdir './ic_auto_selection' 's'
% end

% Now do set up to make future coding easier and faster
% Describe the first (Factor) and second level layers
Factors = struct;
Layer_Factors = {'WholeBrain', 'InOutBrain', 'Networks', 'Regions', 'Power', 'Smoothness', 'Corr', 'Conditions'};
Layer_WholeBrain = {'fullvolume_smoothed_nothresh', 'fullvolume_nothresh', 'fullvolume_noclustering', 'fullvolume'};
Layer_InOutBrain = {'Inbrain', 'Outbrain'};
Layer_Networks = {'MedialVisual', 'SensoryMotor', 'DorsalAttention', 'VentralAttention', 'FrontoParietal', 'DefaultModeNetwork', 'Subcortical'};
Layer_Regions = {'Signal', 'GM', 'WM', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'}; % all smaller regions
Layer_Power = {'Lowfreq', 'BOLDfreq', 'Highfreq'};
Layer_Smoothness = {'Smoothing_Retention', 'Clustering_Prop', 'Smoothness'};
Layer_Corr = {'DVARS', 'FD'};
if ~strcmp(task_events_file, '')
    Layer_Conditions = {'General', 'Task'}; % if we have conditions, use them
else
    Layer_Conditions = {'General'}; % otherwise you only have general
end

Layers_Combined = {Layer_WholeBrain, Layer_InOutBrain, Layer_Networks, Layer_Regions, ...
    Layer_Power, Layer_Smoothness, Layer_Corr, Layer_Conditions};

% Now Make struct from the Layers
for i1 = 1:length(Layers_Combined)
    c = cell(length(Layers_Combined{i1}),1);
    Factors.(Layer_Factors{i1}) = cell2struct(c, Layers_Combined{i1});
end

% Also create ready made tags for easier future selection: For First Layer
Tag_Vals = {'WholeBrain', 'InOutBrain', 'Networks', 'Regions'}; % which ones you will calculate ICsum for, in Factor Layer
Tag_Prop_Regions = {'Networks', 'Regions'}; % Which ones you will calculate proportions for (by dividing by regions_compare)
Tag_Prop_InOut = {'InOut'}; % Which ones you will calculate proportions for (by dividing by InOut Total)
Tag_Prop_Power = {'Power'}; % Which ones you will calculate proportions for (by dividing by total power)
Tag_kmeans = {'InOutBrain', 'Networks', 'Regions', 'Power', 'Smoothness', 'Corr', 'Conditions'}; % which ones we calculate kmeans for
% keep in mind you divide by signal prop if it is Regions, and then
% condition corr is one level deeper than the others

% For Second/Other Layers
Tag_Regions_Compare = {'Signal', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'}; % which regions will we compare all against
Tag_Good = {'Inbrain', 'Signal', 'GM', 'BOLDfreq', 'HRF', 'Smoothness'}; % High kmeans is signal
Tag_Bad = {'Outbrain', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly', 'Lowfreq', 'Highfreq', 'DVARS_corr', 'Move_corr'}; % High kmeans is noise


% (1) Check Spatial Map Overlap with ROIs
% (1a) Approximate overlap with mean and numvoxels

% We can just make tables for everything too
WholeBrain_labels = {'fullvolume_smoothed_nothresh', 'fullvolume_nothresh', 'fullvolume_noclustering', 'fullvolume'};
ROI_labels = {'Signal', 'GM', 'WM', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'};
InOut_labels = {'Inbrain', 'Outbrain'};
Freq_labels = {'Lowfreq', 'BOLDfreq', 'Highfreq'};
Corr_labels = {'FD_Corr', 'DVARS_Corr'};
Smoothness_labels = {'Smoothing_Retention', 'Clustering_Prop', 'Smoothness'};
Networks_labels = {'MedialVisual', 'SensoryMotor', 'DorsalAttention', 'VentralAttention', 'FrontoParietal', 'DefaultModeNetwork', 'Subcortical'};
if isfile(task_events_file)
    HRF_labels = {'HRF_general', 'HRF_task'};
else
    HRF_labels = {'HRF_general'};
end

% Calculate ICsums for relevant things
ICsum_labels = [WholeBrain_labels, InOut_labels, Networks_labels, ROI_labels];
ICmean_table = table();
ICnumvoxels_table = table();
ICsum_table = table();
for i1 = 1:length(ICsum_labels)
    ICmean_table{:,ICsum_labels{i1}} = dlmread(['ROIcalcs/', ICsum_labels{i1}, '_ICmean.txt']);
    ICnumvoxels_table{:,ICsum_labels{i1}} = dlmread(['ROIcalcs/', ICsum_labels{i1}, '_ICnumvoxels.txt']);
    ICsum_table{:,ICsum_labels{i1}} =  ICmean_table{:,ICsum_labels{i1}} .* ICnumvoxels_table{:,ICsum_labels{i1}};
end
Tables.ICmean = ICmean_table; Tables.ICnumvoxels = ICnumvoxels_table; Tables.ICsum_table = ICsum_table;

% Now Calculate the proportaions for relevant things:
% For ROIs and Networks, just divide by the following:
ROI_compare_labels = {'Signal', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'};
ROINetworks_labels = [ROI_labels, Networks_labels];
ROINetworks_genprops_table = table();
for i1 = 1:length(ROINetworks_labels)
    ROINetworks_genprops_table{:,ROINetworks_labels{i1}} = ICsum_table{:,ROINetworks_labels{i1}} ./ sum(ICsum_table{:, ROI_compare_labels},2);
end
ROI_genprops_table = table(); ROI_genprops_table{:, ROI_labels} = ROINetworks_genprops_table{:,ROI_labels};
Networks_genprops_table = table(); Networks_genprops_table{:, Networks_labels} = ROINetworks_genprops_table{:,Networks_labels};
Tables.ROI_generalprops_table = ROI_genprops_table;
Tables.Networks_generalprops_table = Networks_genprops_table;

% For ROI Noise areas, we will make decisions in comparison to signal
% overlap for the IC
ROI_relprops_labels = {'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'}; % Noise parts only
ROI_relprops_table = table();
ROI_relprops_table{:,ROI_labels} = ROINetworks_genprops_table{:,ROI_labels};
% Recalculate noise compartments by dividing by signal:
for i1 = 1:length(ROI_relprops_labels)
    ROI_relprops_table{:, ROI_relprops_labels{i1}} = ROINetworks_genprops_table{:,ROI_relprops_labels{i1}} ...
        ./ (ROINetworks_genprops_table{:,'Signal'} + ROINetworks_genprops_table{:,ROI_relprops_labels{i1}});
end
Tables.ROI_relativeprops_table = ROI_relprops_table;

% Inbrain Outbrain should just be compared to themselves for gen props
InOut_compare_labels = {'Inbrain', 'Outbrain'};
InOut_genprops_table = table();
for i1 = 1:length(InOut_labels)
    InOut_genprops_table{:,InOut_labels{i1}} = ICsum_table{:,InOut_labels{i1}} ./ sum(ICsum_table{:, InOut_compare_labels},2);
end
Tables.InOutBrain_generalprops_table = InOut_genprops_table;

% Outbrain relative props can be based on signal as well
InOut_relprops_table = InOut_genprops_table;
InOut_relprops_table{:, 'Outbrain'} = InOut_relprops_table{:, 'Outbrain'} ./ (ROINetworks_genprops_table{:,'Signal'} + InOut_relprops_table{:, 'Outbrain'});
Tables.InOut_relativeprops_table = InOut_relprops_table;

% While we are here, let's initialize IC number information (we needed to
% dlmread some files first)
% Grab number of ICs and explained variance:
ICs = [1:length(ROI_genprops_table{:,'Signal'})]';
ICs_table = array2table(ICs, 'VariableNames', {'ICs'});
IC_exp_var = dlmread('ROIcalcs/IC_exp_variance.txt');
IC_exp_var = IC_exp_var ./ 100; % put it in decimal places
IC_exp_var_table = array2table(IC_exp_var, 'VariableNames', {'Variance_Explained'});
% Add to Tables Struct
Tables.IC_table = ICs_table;
Tables.IC_exp_var_table = IC_exp_var_table;

% And also calculate smoothing values
Smoothness_table = table();
Smoothness_table{:, 'Smoothing_Retention'} = ICsum_table{:,'fullvolume_smoothed_nothresh'} ./ ICsum_table{:, 'fullvolume_nothresh'}; % this is the one we use
Smoothness_table{:, 'Clustering_Prop'} = ICsum_table{:,'fullvolume'} ./ ICsum_table{:, 'fullvolume_noclustering'};
Smoothness_table{:, 'Smoothness'} = 1 ./ dlmread('clustering/smoothness.txt');
Tables.Smoothness_table = Smoothness_table;

% (2) Check Power Frequency Analysis in low, BOLD (0.008-0.15), and high
% frequencies
% (2a) Grab time series and calculate power and proportions in each
% frequency range of interest
N = T/dt; F = 1/dt;
df = 1/T; N_freq = N/2 + 1;
f = F*(0:floor(N/2))/N; % to chart what frequencies we are at
frequencies = f; % easier naming later
lower_phys_cutoff = round(0.008 / df) + 1; % start at freq 0 at position 1
higher_phys_cutoff = round(0.15 / df) + 1;
cush = 1; % to avoid overlap, give an index of padding

% while we are here, let's also calculate hrf general powerspectrum
Data.HRF_general.plot_norm = normalize(hrf_padded)'; % remove mean to get rid of 0Hz
Data.HRF_general.fft = fft(Data.HRF_general.plot_norm);
Data.HRF_general.P2_single_hrf = abs((Data.HRF_general.fft.^2)/N); % power
Data.HRF_general.hrf_power = Data.HRF_general.P2_single_hrf(1:floor(N/2)+1); % grab first half of power, since power has negative to positive frequency
% normalize so it sums to 1
Data.HRF_general.hrf_power_norm = (Data.HRF_general.hrf_power ./ trapz(Data.HRF_general.hrf_power)); % if you multiply this by the normalized power for each IC, you may be able to calculate useful things
P1_single_hrf_bp = Data.HRF_general.hrf_power_norm;
P1_single_hrf_bp(1:lower_phys_cutoff) = 0;
P1_single_hrf_bp(higher_phys_cutoff+1:length(f)) = 0; % bp is bandpassed
Data.HRF_general.hrf_power_bp = P1_single_hrf_bp;
% Can plot(HRF_general.f,HRF_general.hrf_power_bp) to see the power if you want to confirm

% OK, now calculate power spectrums for each IC:
% set size of arrays explicitely to save computation time
ts = zeros(numvolumes, length(ICs));
Ps = zeros(numvolumes/2+1, length(ICs)); % Powerspectrums
f_all_power = zeros(length(ICs), 1);
BOLDfreqIC = zeros(length(ICs), 1);
LowfreqIC = zeros(length(ICs), 1);
HighfreqIC = zeros(length(ICs), 1);
for i=1:length(ICs)
    % grab time spectrum
    file=strcat('melodic/report/t', num2str(i), '.txt');
    ts(:,i)=dlmread(file);

    % calculate power spectrum as evidenced in matlab tutorials (floor is
    % used to handle both even and odd input lengths of N)
    t = ts(:,i);
    Y = fft(t);
    P2 = abs((Y.^2)/N); % power
    P1 = P2(1:floor(N/2)+1); % grab first half of power, since power has negative to positive frequency
    f = F*(0:floor(N/2))/N; % to chart what frequencies we are at
    P1_norm = P1 ./ trapz(P1);
    Ps(:,i) = P1_norm; % so it sums to 1 for easier comparison to hrf

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
% Add some of this information into a Data Struct for later reference
Data.Powers = Ps;
Data.TimeSeries = ts;
Data.N = N; Data.T = T; Data.tr = dt; Data.frequencies = f;
Data.lowfreq_cutoff = 0.008; Data.highfreq_cutoff = 0.15;

% Calculate overlap similarity and record in table, this may be added to
% later
HRF_table = table();
HRF_table{:,'general_power_overlap'} = normalize(sum(Data.HRF_general.hrf_power_bp .* Data.Powers), "range")'; % a higher sum for an IC indicates more power spectrum overlap with expected BOLD response
Tables.HRF_table = HRF_table;

% (2b) Calculate actual measured proportions, easy to do manually.
LowfreqIC_prop = LowfreqIC ./ (LowfreqIC+BOLDfreqIC+HighfreqIC);
BOLDfreqIC_prop = BOLDfreqIC ./ (LowfreqIC+BOLDfreqIC+HighfreqIC);
HighfreqIC_prop = HighfreqIC ./ (LowfreqIC+BOLDfreqIC+HighfreqIC);

% Add to a frequency table
Freq_genprops_table = table();
Freq_genprops_table{:,Freq_labels} = [LowfreqIC_prop, BOLDfreqIC_prop, HighfreqIC_prop];
Tables.Freq_generalprops_table = Freq_genprops_table;

% And make a relprops version where you divide by BOLDfreqIC_prop, can do
% without a loop since it is so short
Freq_relprops_table = Freq_genprops_table;
Freq_relprops_table{:,'Lowfreq'} = Freq_genprops_table{:,'Lowfreq'} ./ (Freq_genprops_table{:, 'BOLDfreq'} + Freq_genprops_table{:,'Lowfreq'});
Freq_relprops_table{:,'Highfreq'} = Freq_genprops_table{:,'Highfreq'} ./ (Freq_genprops_table{:, 'BOLDfreq'} + Freq_genprops_table{:,'Highfreq'});
Tables.Freq_relativeprops_table = Freq_relprops_table;

% (3) Correlations: to dvars, FD, tasks (if not resting state):
confound_place = 'confounds_timeseries.csv';
allconfounds = readtable(confound_place);
Data.Confounds.Nineparam = table2array(allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'white_matter', 'csf', 'global_signal'}));
Data.Confounds.Eightparam = table2array(allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z', 'white_matter', 'csf'}));
Data.Confounds.DVARS = table2array(allconfounds(:,{'dvars'}));
Data.Confounds.FD = table2array(allconfounds(:,{'framewise_displacement'}));

% Calculate and create Correlations table
Corr_genprops_table = table();
% FD: 
corr_during_FD = corr(detrend(Data.Confounds.FD(2:end)), abs(diff(detrend(Data.TimeSeries)))).^2';
Corr_genprops_table{:, 'FD_Corr'} = corr_during_FD;
% DVARS: 
corr_during_DVARS = corr(detrend(Data.Confounds.DVARS(2:end)), abs(diff(detrend(Data.TimeSeries)))).^2';
Corr_genprops_table{:, 'DVARS_Corr'} =  corr_during_DVARS;
Tables.Corr_generalprops_table = Corr_genprops_table;


% calculate condition hrfs if a task_events_file was provided, can be
% helpful in deciding if an IC is task related
condition_names = {};
condition_corrs = [];
if isfile(task_events_file)
    condition_names = unique(task_events_table{:,3})';
    block = zeros(numvolumes, length(condition_names));
    
    for j = 1:length(condition_names)
        curr_condition_name = condition_names{j};
        curr_condition_rows = find(strcmp(task_events_table{:,3}, curr_condition_name));
        for k = 1:length(curr_condition_rows)
            start_ind = task_events_table{curr_condition_rows(k),1}./TR+1;
            stop_ind = start_ind + task_events_table{curr_condition_rows(k),2}./TR - 1;
            block(start_ind:stop_ind, j) = 1;
        end
    end

    events_table = array2table(block, 'VariableNames', condition_names);

    % The above should give us blocks for all conditions, now convolve with hrf
    % to get a guestimated hrf response for each condition:
    hrf_conditions_table = table();
    hrf_conditions = zeros(size(block));
    for j = 1:length(condition_names)
        convolution = conv(Data.HRF_general.plot_norm, block(:,j));
        hrf_conditions(:, j) = convolution(1:size(block,1));
    end
   
    condition_names_orig = condition_names; % save original
    hrf_conditions(:,length(condition_names_orig)+1) = sum(hrf_conditions,2); 
    hrf_conditions = normalize(hrf_conditions); % normalize to remove mean so we don't have power at 0Hz overtaking
    condition_names = [condition_names, 'Combined']; 
    condition_corrs = zeros(length(ICs), length(condition_names));

    condition_corrs_table = table();
    for m = 1:length(condition_names)
        condition_corrs(:,m) = corr(hrf_conditions(:,m), ts).^2;
    end

    % make hrf_conditions into a table
    hrf_conditions_table = array2table(hrf_conditions, 'VariableNames', condition_names);
    Tables.hrf_conditions_table = hrf_conditions_table;

    % make condition_corrs table
    condition_corrs_table = array2table(condition_corrs, 'VariableNames', condition_names);
    Tables.condition_corrs_table = condition_corrs_table;

    % Save into Data struct as appropriate
    Data.Task.EventBlocks = events_table;
    Data.Task.hrf_conditions = hrf_conditions_table;
end



% And then what about for task correlated HRF?
% condition correlated, if we have it (and also make a power spectrum of
% most correlated condition corr)
if ~isempty(condition_corrs)
    % repeat IC matrix to match size
    ICs_reps = repmat(ICs, 1, size(condition_corrs,2));
    
    % both calculate relative correlation strength for hrfs AND calculate
    % power spectra
    % initialize things first, then go into for loop
    N = T/dt; F = 1/dt; df = 1/T; N_freq = N/2 + 1;
    Y_all_hrfs = zeros(length(ts), length(condition_names));
    P2_all_hrfs = Y_all_hrfs;
    P1_all_hrfs = zeros(length(1:floor(N/2)+1), length(condition_names));
    P1_all_hrfs_norm = P1_all_hrfs; P1_all_hrfs_bp = P1_all_hrfs_norm;
    for m = 1:length(condition_names)
        % power spectra
        curr_hrf = hrf_conditions(:,m);
        Y_all_hrfs(:,m) = fft(curr_hrf);
        P2_all_hrfs(:,m) = abs((Y_all_hrfs(:,m).^2)/N); % power
        P1_all_hrfs(:,m) = P2_all_hrfs(1:floor(N/2)+1,m); % grab first half of power, since power has negative to positive frequency
        % normalize so it sums to 1
        P1_all_hrfs_norm(:,m) = (P1_all_hrfs(:,m) ./ trapz(P1_all_hrfs(:,m))); % if you multiply this by the normalized power for each IC, you may be able to calculate useful things
        P1_all_hrfs_bp(:,m) = P1_all_hrfs_norm(:,m);
        P1_all_hrfs_bp(1:lower_phys_cutoff,m) = 0;
        P1_all_hrfs_bp(higher_phys_cutoff+1:length(f),m) = 0; % bp is bandpassed
        % Can plot(f,P1_all_hrfs_norm) to see the power if you want
    end

    % Add to relevant struct:
    Data.HRF_task.plot_norm = hrf_conditions; % Was normalized earlier
    Data.HRF_task.fft = Y_all_hrfs;
    Data.HRF_task.P2_single_hrf = P2_all_hrfs;
    Data.HRF_task.hrf_power = P1_all_hrfs;
    Data.HRF_task.hrf_power_norm = P1_all_hrfs_norm;
    Data.HRF_task.hrf_power_bp = P1_all_hrfs_bp;

    % Use the maximum correlation to pick out the best hrf_power spectrum
    % to use per IC
    [M,I] = max(condition_corrs, [], 2);
    condition_corr_most_correlated_index = I;
    max_conditioncorr_powerspectra = P1_all_hrfs_bp(:,I)'; % each row is an IC, each column is a positive powerspectrum

    % OK, NOW we can calculate relative poweroverlap
    power_overlap_all_hrfs = normalize(sum(max_conditioncorr_powerspectra .* Ps', 2), "range"); % a higher sum for an IC indicates more power spectrum overlap with expected BOLD response
    hrfstask_powerspectrum = P1_all_hrfs_norm;
    hrfstask_fullpowerspectrum = P2_all_hrfs;
    HRF_table{:, 'besttask_power_overlap'} = power_overlap_all_hrfs;

    
    save('Task_Condition_Info.mat', 'ts', 'numvolumes','hrf_conditions_table', 'condition_corr_most_correlated_index', ...
        'condition_corrs', 'condition_names', 'hrfstask_fullpowerspectrum', ...
        'hrfstask_powerspectrum', 'max_conditioncorr_powerspectra', 'power_overlap_all_hrfs', 'frequencies', 'hrf_conditions')
end

% save standard 9Param and 8Param for comparison or use later:
writematrix(Data.Confounds.Nineparam, ['9p_regressors.txt'], 'Delimiter', ' ')
writematrix(Data.Confounds.Eightparam, ['8p_regressors.txt'], 'Delimiter', ' ')

% Calculate the highest normalized HRF matching (whether to a general HRF
% frequency range, or HRFs in response to tasks)
HRF_table{:, 'best_power_overlap'} = max(HRF_table{:,:}, [], 2); % if there is no task, it will just be the same as general
Tables.HRF_table = HRF_table;


% (4) Putting it all together!

% Make a normalized array so it is easier to follow in the future

% the general values are true proportions
feature_general_table = [ICs_table, ROI_genprops_table, InOut_genprops_table, ...
    Freq_genprops_table, Corr_genprops_table, Smoothness_table, HRF_table, ...
    Networks_genprops_table, condition_corrs_table];
% relative table is the one we will use for kmeans things:
feature_relative_table = [ICs_table, ROI_relprops_table, InOut_relprops_table, ...
    Freq_relprops_table, Corr_genprops_table, Smoothness_table, HRF_table, ...
    Networks_genprops_table, condition_corrs_table];
% norm table can put things in perspective for use
feature_table_norm = [feature_relative_table(:,1), normalize(feature_relative_table(:,2:end))]; % normalized version
feature_labels = feature_relative_table.Properties.VariableNames;

% Add to Tables
Tables.feature_general_table = feature_general_table;
Tables.feature_relative_table = feature_relative_table;
Tables.feature_table_norm = feature_table_norm;


% Label everything based on kmeans! Nice and easy with tables
classification_table = table();
classification_table{:, 'ICs'} = ICs;
for i1 = 2:length(feature_labels)
    curr_val = feature_relative_table{:,feature_labels{i1}};
    [val_idx, val_C] = kmeans(curr_val, 3, 'Start', [min(curr_val); median(curr_val); max(curr_val)]);
    classification_table{:, ['High_', feature_labels{i1}]} = val_idx == find(val_C == max(val_C));
    classification_table{:, ['Low_', feature_labels{i1}]} = val_idx == find(val_C == min(val_C));
end

% Start parsing out what is important to look at from classification
best_classification_labels = {'High_Signal', 'High_Smoothing_Retention', 'High_best_power_overlap'};
othergood_classification_labels = {'High_GM', 'High_BOLDfreq', 'High_general_power_overlap', 'High_besttask_power_overlap', 'High_WM'};
bad_classification_labels = {'Low_Signal', 'Low_GM', 'High_Edge', 'High_Subepe', 'High_CSF', 'High_Suscept', ...
    'High_OutbrainOnly', 'High_Outbrain', 'High_Lowfreq', 'Low_BOLDfreq', 'High_Highfreq', ...
    'Low_best_power_overlap', 'Low_Smoothing_Retention', 'High_DVARS_Corr', 'High_FD_Corr'}; 
networks_classification_labels = {'High_MedialVisual', 'High_SensoryMotor', ...
    'High_DorsalAttention', 'High_VentralAttention', 'High_FrontoParietal', ...
    'High_DefaultModeNetwork', 'High_Subcortical'};

for i1 = 1:length(condition_names)
    condition_corrs_classification_labels{i1} = ['High_', condition_names{i1}];
end

best_classification_table = classification_table(:, best_classification_labels);
othergood_classification_table = classification_table(:, othergood_classification_labels);
good_classification_table = [best_classification_table, othergood_classification_table];
bad_classification_table = classification_table(:, bad_classification_labels);
networks_classification_table = classification_table(:, networks_classification_labels);
condition_corrs_classification_table = classification_table(:, condition_corrs_classification_labels);

% remake classification in a condensed format
classification_final_table = [ICs_table, good_classification_table, bad_classification_table, networks_classification_table, condition_corrs_classification_table];
classification_general_table = [ICs_table, good_classification_table, bad_classification_table]; % good for general labeling

% Add all to tables:
Tables.classification_table = classification_table;
Tables.classification_final_table = classification_final_table;
Tables.classification_general_table = classification_general_table;
Tables.best_classification_table = best_classification_table;
Tables.othergood_classification_table = othergood_classification_table;
Tables.good_classification_table = good_classification_table;
Tables.bad_classification_table = bad_classification_table;
Tables.networks_classification_table = networks_classification_table;
Tables.condition_corrs_classification_table = condition_corrs_classification_table;

% Take into account low and high frequency
% content into account for signal and noise
% high signal
% intialize Results Struct
Results = struct();
% high noise
Highnoise_indices = sum(table2array(bad_classification_table), 2) > 0;
Highnoise_ICs = ICs(Highnoise_indices);

% all 3 of the best things are classified high
HighSig_indices = (sum(table2array(best_classification_table), 2) == width(best_classification_table)); 
HighSig_ICs = ICs(HighSig_indices);

% Also make something to mark whether it has any good classification tags
HasSig_indices = sum(table2array(good_classification_table),2) > 0;
HasSig_ICs = ICs(HasSig_indices);

% Add these to results struct
Results.Highnoise_indices = Highnoise_indices;
Results.Highnoise_ICs = Highnoise_ICs;
Results.HasSig_indices = HasSig_indices;
Results.HasSig_ICs = HasSig_ICs;
Results.HighSig_indices = HighSig_indices;
Results.HighSig_ICs = HighSig_ICs;

% OK, now we can start to narrow down our search! Good ICs are
% determined by GM overlap (compared to other regions),
% Boldfrequency power (compared to high frequency), and
% smoothness. Weight GM overlap heavier than smoothness or
% frequency (because noise could be in BOLD and smooth, like
% WMCSF signal)
[signal_val, signal_idx] = sort(normalize(feature_relative_table{:,'Smoothing_Retention'}, "range") .* ...
    normalize(feature_relative_table{:,'Signal'}, "range").^2 .* ...
    normalize(feature_relative_table{:,'best_power_overlap'}, "range"), 'descend');

% A non-smooth IC can be saved with high network overlap if HRF and
% signal/GM prop values are also relatively high. This can help save us sometimes
% (e.g., subcortical is often not very smooth, so could use saving)
% Not all networks may be found, so it is likely best to do a new kmeans
% that looks at all network overlap:
network_props = feature_relative_table{:, Layer_Networks};
[rows, cols] = size(network_props);
curr_val = reshape(network_props, [rows.*cols, 1]);
[val_idx, val_C] = kmeans(curr_val, 3, 'Start', [min(curr_val); median(curr_val); max(curr_val)]);
high_network_props = val_idx == find(val_C == max(val_C));
high_network_props = reshape(high_network_props, size(network_props));

Network_Overlap_indices = (sum(high_network_props, 2) > 0) & ...
    (sum([classification_final_table{:, 'High_Signal'}, classification_final_table{:, 'High_GM'}],2) > 0) ...
    & (classification_final_table{:, 'High_best_power_overlap'} == 1);
Network_Overlap_ICs = ICs(Network_Overlap_indices);

% Add to results struct
Results.Network_Overlap_indices = Network_Overlap_indices;
Results.Network_Overlap_ICs = Network_Overlap_ICs;
% still display all the high networks, just for reference if ever
% needed/desired

% Each time you get noise, drop one in tolerance, if you get signal, add
% one in tolerance (up to tolerance number), stop when tolerance is 0!
g = 1;
signal_decision = [];
stop_num = tolerance;
while stop_num > 0
    % if it is labeled as high noise, it could be saved if it has good
    % Network Overlap and great HRF overlap
    if ismember(signal_idx(g), Highnoise_ICs) == 1
        if ismember(signal_idx(g), Network_Overlap_ICs) == 0
            stop_num = stop_num - 1;
            % mark as noise!
            signal_decision(g) = 0;
        else
            stop_num = stop_num + 1;
            % mark as signal
            signal_decision(g) = 1;
        end
    elseif ismember(signal_idx(g), HasSig_ICs) == 1
        % Needs to have at least one signal-like marking
        stop_num = stop_num + 1;
        % mark as signal
        signal_decision(g) = 1;
    else
        % if it has no signal-like markings, mark as noise
        stop_num = stop_num - 1;
        % mark as signal
        signal_decision(g) = 0;
    end
       % if stop num went above tolerance, set it back to tolerance
       if stop_num > tolerance
           stop_num = tolerance;
       end
    g = g + 1; % g will be how far along the signal ICs you go
end
g = g - 1;

% Make relevant table and give appropriate signal labeling:
ICs_signal_chance_order_table = array2table(signal_idx, 'VariableNames', {'ICs'});

potential_signal_ICs = signal_idx(1:g);
signal_indices = zeros(1,length(ICs));
signal_indices(potential_signal_ICs) = signal_decision; % Labels based on signal decision earlier, and anything past that as noise
signal_indices = logical(signal_indices');


% Give a fail safe incase it is so bad that no ICs are labeled as signal
if sum(signal_indices) == 0
    fprintf('WARNING: No ICs were labeled as signal. Will label first IC in best order as signal to allow the rest of the code to run.\n')
    signal_indices(:,1) = logical(1);
end
noise_indices = ~signal_indices;
signal_ICs = ICs(signal_indices);
noise_ICs = ICs(noise_indices);

potential_signal_ICs = signal_idx; % Just make sure that you include all ICs all so that we can just always show all of them
signal_checker_table = array2table(potential_signal_ICs, 'VariableNames', {'PotentialICs'});
signal_checker_table.SignalLabel = signal_indices(potential_signal_ICs);
signal_checker_table.HighSignalLabel = HighSig_indices(potential_signal_ICs);
signal_checker_table.HighNoiseLabel = Highnoise_indices(potential_signal_ICs);
Highconditioncorr_indices = sum(table2array(condition_corrs_classification_table),2) > 0;
if ~isempty(condition_names)
    signal_checker_table.TaskCorrelated = double(sum(Highconditioncorr_indices(potential_signal_ICs,:),2) > 0);
end


% add on classification labels to signal_checker_table 
% One column for good, one for bad, one for networks, one for tasks (if it
% exists)
% good first:
good_classification_potentialsignal_labels = cell(length(potential_signal_ICs), 1);
for n = 1:length(potential_signal_ICs)
    curr_class = {good_classification_table.Properties.VariableNames{logical(good_classification_table{potential_signal_ICs(n),:})}};
    % put 'None' if none exist, helps visualize on excel
    if isempty(curr_class)
        curr_class = {'None'};
    end
    good_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
end
signal_checker_table.Good_Tags = good_classification_potentialsignal_labels;

% now bad:
bad_classification_potentialsignal_labels = cell(length(potential_signal_ICs), 1);
for n = 1:length(potential_signal_ICs)
    curr_class = {bad_classification_table.Properties.VariableNames{logical(bad_classification_table{potential_signal_ICs(n),:})}};
    % put 'None' if none exist, helps visualize on excel
    if isempty(curr_class)
        curr_class = {'None'};
    end
    bad_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
end
signal_checker_table.Bad_Tags = bad_classification_potentialsignal_labels;

% and networks
networks_classification_potentialsignal_labels = cell(length(potential_signal_ICs), 1);
for n = 1:length(potential_signal_ICs)
    curr_class = {networks_classification_table.Properties.VariableNames{logical(networks_classification_table{potential_signal_ICs(n),:})}};
    % put 'None' if none exist, helps visualize on excel
    if isempty(curr_class)
        curr_class = {'None'};
    end
    networks_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
end
signal_checker_table.Networks_Tags = networks_classification_potentialsignal_labels;

% finally condition corrs (if they exist)
if ~isempty(condition_names)
    condition_corrs_classification_potentialsignal_labels = cell(length(potential_signal_ICs), 1);
    for n = 1:length(potential_signal_ICs)
        curr_class = {condition_corrs_classification_table.Properties.VariableNames{logical(condition_corrs_classification_table{potential_signal_ICs(n),:})}};
        % put 'None' if none exist, helps visualize on excel
        if isempty(curr_class)
            curr_class = {'None'};
        end
        condition_corrs_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
    end
    signal_checker_table.Condition_Tags = condition_corrs_classification_potentialsignal_labels;
end

IC_checker_table = signal_checker_table; % renamed to something more recognizable

% Add important things to Results struct:
Results.tolerance = tolerance;
Results.IC_checker_table = IC_checker_table;
Results.feature_relative_table = feature_relative_table;
Results.feature_table_norm = feature_table_norm;
Results.signal_indices = signal_indices;
Results.signal_ICs = signal_ICs;
Results.noise_indices = noise_indices;
Results.noise_ICs = noise_ICs;

% If we did manual IC checker through this script, we would update
% IC_checker_table and then run just the part below. But we now have a
% separate script for this. While redundent, this is still good to have as
% a reference if every needed.
%%%%% Update IC_checker_updated_table variable above before
%%%%% re-evaluating everything below
% need to reconvert table to not classification column(s) 
cut_off_table = IC_checker_table(:,{'PotentialICs', 'SignalLabel', 'HighSignalLabel', 'HighNoiseLabel'});
close_IC_checker_updated_array = table2array(cut_off_table);
% update noise_indices
noise_indices(close_IC_checker_updated_array(close_IC_checker_updated_array(:,2) == 0,1)) = 1;
noise_indices(close_IC_checker_updated_array(close_IC_checker_updated_array(:,2) == 1,1)) = 0;
noise_indices = logical(noise_indices);

% Create final noise and signal ICs and indices, if GM_prop>0.5,
% label as signal regardless
% check close ICs and update noise indices
noise_ICs = ICs(noise_indices);
signal_indices = logical(~noise_indices);
signal_ICs = ICs(signal_indices);

% Update Results struct
Results.signal_indices = signal_indices;
Results.signal_ICs = signal_ICs;
Results.noise_indices = noise_indices;
Results.noise_ICs = noise_ICs;

% (3b) Make it easier to compare groups by features
eval_table = feature_table_norm(close_IC_checker_updated_array(:,1), :);
Results.eval_table = eval_table;

% (3d) Compare before and after selection
feature_relative_array = table2array(feature_relative_table(:,2:end));
before = sum(feature_relative_array .* IC_exp_var);
IC_exp_var_signal = IC_exp_var(signal_ICs);
IC_exp_var_signal = IC_exp_var_signal ./ sum(IC_exp_var_signal); % this re-does the proportions for just the signal components
after = sum(feature_relative_array(signal_ICs,:) .* IC_exp_var_signal);
percent_change = 100.*(after ./ before - 1); % rounded to be easier to read
compare_cleaning = array2table([before; after; percent_change]', 'RowNames', feature_relative_table.Properties.VariableNames(2:end), 'VariableNames', {'Before', 'After', 'Percent_Change'});
Results.compare_cleaning = compare_cleaning;

% estimate effective degrees of freedom (somewhat) from CADICA alone via proportions:
Results.percent_variance_kept = sum(IC_exp_var(signal_ICs)); % probably want > 0.1 (10%)
Results.num_ICs_kept = length(signal_ICs);
Results.num_ICs_total = length(ICs);
Reults.IC_exp_var = IC_exp_var;
Results.ICs = ICs;

% Create images to review signal and noise components from
% probabilities ICs
all_prob = niftiread('./melodic/ICprobabilities.nii.gz');
all_clusterthresh = niftiread('./melodic/ICclusterthresh_zstat.nii.gz');

% higher probability clustered & Highest prob only clustered
all_prob_clustered = cast(all_prob .* (all_prob > 0.67) .* (abs(all_clusterthresh) > 3), 'single'); % includes probability values
all_highest_prob_clustered = cast(all_prob .* (all_prob > 0.95) .* (abs(all_clusterthresh) > 3), 'single'); % higher cut offs 

gm_prob = niftiread('./region_masks/GM_prop_final.nii.gz');
notgm_prob = niftiread('./region_masks/NotGM_prop_final.nii.gz');
gm_bin = gm_prob > 0.75;
notgm_bin = notgm_prob > 0.75;
noise_prob = all_prob_clustered(:,:,:,noise_ICs);
noise_prob_info = niftiinfo('./melodic/ICprobabilities.nii.gz');
signal_prob = all_prob_clustered(:,:,:,signal_ICs);

% grab noise probabilities, update header, and write, same for
% signal
noise_prob_info.ImageSize = size(noise_prob);
noise_prob_info.DisplayIntensityRange = [0.9 1];
noise_prob_info.Datatype = 'single';

signal_prob_info = noise_prob_info;
signal_prob_info.ImageSize = size(signal_prob);

potential_signal_prob = all_prob_clustered(:,:,:,potential_signal_ICs);
potential_signal_prob_info = noise_prob_info;
potential_signal_prob_info.ImageSize = size(potential_signal_prob);

niftiwrite(potential_signal_prob, 'PotentialSignalICs', potential_signal_prob_info, 'Compressed', true)
niftiwrite(noise_prob, 'NoiseICs', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob, 'SignalICs', signal_prob_info, 'Compressed', true)

% Now compress noise and signal, with only max probabilities per voxel
[all_prob_1D, all_prob_ind] = max(all_prob, [], 4); 
[noise_prob_1D, noise_prob_ind] = max(all_prob(:,:,:,noise_ICs), [], 4); 
[signal_prob_1D, signal_prob_ind] = max(all_prob(:,:,:,signal_ICs), [], 4); 

% Update headers for the 3D images
noise_prob_info.ImageSize = noise_prob_info.ImageSize(1:end-1);
noise_prob_info.PixelDimensions = noise_prob_info.PixelDimensions(1:end-1);
signal_prob_info = noise_prob_info;

% Calculate out a signal to noise ratio mask overlap - but this might not
% be very helpful
signal_noise_ratio_IC_overlap = signal_prob_1D ./ (noise_prob_1D+signal_prob_1D);
signal_noise_ratio_IC_overlap(isnan(signal_noise_ratio_IC_overlap)) = 0.5; % convert if/when there is no noise or signal high prob value to 0.5 (will become 0)
signal_noise_ratio_IC_overlap = signal_noise_ratio_IC_overlap - 0.5; % shift 0.5 to 0, to center it so <0 is more nois & >0 is more signal

% perhaps a better method is threshold at 0.95 both signal and noise, make
% noise be noise - signal (to get noise sections alone), set to -1 for
% noise only sections, then combine signal and noise masks together.
thresholded_noise = noise_prob_1D > 0.949;
thresholded_signal = signal_prob_1D > 0.949;
noise_alone = thresholded_noise .* ~thresholded_signal;
signal_and_noise_overlap = thresholded_signal - noise_alone; % should just be 1 for signal, and -1 for noise

% write out the masks of highest probability
niftiwrite(noise_prob_1D, 'NoiseICOverlap', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob_1D, 'SignalICOverlap', signal_prob_info, 'Compressed', true)
niftiwrite(signal_noise_ratio_IC_overlap, 'SignaltoNoiseICOverlap', signal_prob_info, 'Compressed', true) % <0 is more noise, >0 is more signal, 0 is either equivalent, or not high probability either way
niftiwrite(cast(signal_noise_ratio_IC_overlap .* gm_bin, 'single'), 'SignaltoNoiseICOverlap_GM', signal_prob_info, 'Compressed', true) % <0 is more noise, >0 is more signal
niftiwrite(cast(signal_and_noise_overlap, 'single'), 'SignalandNoiseICOverlap', signal_prob_info, 'Compressed', true) % This is likely the most helpful one

% Write out the max noise ICs where we are in gray matter and noise is
% currently more represented, and vice versa. Can help in identifying
% misrepresnted signal and/or noise
int_prob_info = noise_prob_info;
int_prob_info.Datatype = 'int8';
niftiwrite(cast(noise_ICs(noise_prob_ind) .* gm_bin .* (signal_noise_ratio_IC_overlap < 0), 'int8'), 'Max_NoiseIC_GM', int_prob_info, 'Compressed', true)
niftiwrite(cast(signal_ICs(signal_prob_ind).*notgm_bin.*(signal_noise_ratio_IC_overlap>0), 'int8'), 'Max_SignalIC_NotGM', int_prob_info, 'Compressed', true)

% (5c) Export variables later use in next steps (e.g., fsl_regfilt)
% Auto prefix tells you that it is using the automatic best
% guess, and not with manual checking. This should be good, but
% not perfect.
writematrix(Results.noise_ICs', 'auto_noise_dist_ICs.csv')
writematrix(Results.signal_ICs', 'auto_signal_dist_ICs.csv')
writematrix(Results.noise_indices, 'auto_noise_dist_indices.csv')
writematrix(Results.signal_indices, 'auto_signal_dist_indices.csv')

% and also the IC Checker table, normalized feature table, and compare
% cleaning estimates
writetable(IC_checker_table, 'IC_auto_checker.csv')
writetable(feature_relative_table, 'feature_auto_vals.csv')
writetable(feature_table_norm, 'feature_auto_norms.csv')
writetable(compare_cleaning, 'compare_auto_cleaning.csv')

% and then make noise components if you want to do aggressive denoising like in CONN (regression)
mixing_matrix = dlmread('./melodic/melodic_mix');
auto_noise_dist_covariates = mixing_matrix(:, Results.noise_ICs');
writematrix(auto_noise_dist_covariates, 'auto_noise_dist_covariates.txt', 'Delimiter', ' ')
writematrix(auto_noise_dist_covariates, 'auto_noise_dist_covariates.csv')

movefile auto* ic_auto_selection
movefile *SignalIC* ic_auto_selection
movefile *NoiseIC* ic_auto_selection
movefile *IC_auto_checker.* ic_auto_selection
movefile *_auto_* ic_auto_selection
movefile *regressors.* regressors_timeseries
movefile *timeseries.* regressors_timeseries
% move confounds back out though so it is easier to find
cd regressors_timeseries
movefile confounds_timeseries* ../
cd ../

% save a matrix of relevant variables if you want to examine later!
save('DecisionVariables_Auto.mat', 'cadicafol', 'subject_id', 'session_id', 'taskname', 'ICs', 'Data', 'Tables', 'Results')

movefile *Auto* ic_auto_selection

end