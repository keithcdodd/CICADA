function [cleaned_file] = CICADA_2_AutoLabeling(output_dir, task_events_file, tolerance)
% Run this after running CICADA_1 script
% output_dir must be the same one as for the first script, so that this can
% find all the relevant files. Again, it is best if this is in the format
% of CICADA_dir/subj-id/ses-id/task-id (so, a task_dir, really)

% task_events_file should be given if evaluating a task -- will increase IC
% classification accuracy. Task Events file should be the same structure as
% standard bids or fMRIPrep format .csv 

% standard tolerance of 4 is a good baseline. If you increase that number,
% you may label more ICs as signal. Decreasing it may label less ICs as signal.

% General Method of This Script:s
% This script evaluates each IC to determine if it is more likely to be signal or noise. 
% ICs are first ordered from most to least likely to be signal. This is
% determined based on GM spatial overlap, powerspectrum overlap with hrf
% predicted power, and spatial smoothness. K means clustering helps
% identify different types of noise, following manual guidelines.

% easiest way to make sure fsl will work is to make sure call_fsl.m is on
% your path! This is given in fsl installations: e.g., /usr/local/fsl/etc/matlab/call_fsl.m
% You may need to add to startup.m the calls that fsl suggests in their
% documentation when you search "running fsl in matlab" online -
% FslMatlabConfiguration

% check output_dir exists 
if ~isfolder(output_dir)
    fprintf('ERROR: output_dir does not exist!')
    return;
end

% catch if you say there is no task_events_file with 'x'
% check optional inputs
if (exist('task_events_file', 'var') && strcmp(task_events_file, 'x')) || (~exist('task_events_file', 'var')) || ~ischar(task_events_file)
    task_events_file = '';
    fprintf('   Not using a task events file! \n')
else
    fprintf(['  Using task events file: ', task_events_file, '\n'])
end
if ~exist('tolerance', 'var')
    tolerance = 4;
end
fprintf(['   Tolerance is ' num2str(tolerance), '. Standard is 4.\n'])

% read in the functional file to get necessary information
cd(output_dir)
% read in task name, session name, and subject name based on default folder
% structure
[~,task_id,~]=fileparts(pwd);
cd ../
[~,session_id,~]=fileparts(pwd);
cd ../
[~,subject_id,~]=fileparts(pwd);
cd(output_dir)
funcfilename = 'funcfile';
funcfile_info = niftiinfo([funcfilename, '.nii.gz']);

% initialize
Tables = struct;
Data = struct;
Data.tr = funcfile_info.PixelDimensions(4);
Data.numvolumes = funcfile_info.ImageSize(4);
repetitiontime = Data.tr;
numvolumes = Data.numvolumes;

T = repetitiontime .* numvolumes;
TR = repetitiontime;

hrf_conditions_table = table([]);


% model a double gamma hrf, based on spm (but without the need to install
% spm). Still, cite spm since it is completely based on their code.
p = [6 16 1 1 6 0 32]; fMRI_T = 16;
dt = TR/fMRI_T; u = (0:(p(7)/dt)) - p(6)/dt;
h1=p(1)/p(3); l1=dt/p(3);
h2=p(2)/p(4); l2=dt/p(4);
f1 = exp( (h1-1).*log(u) +h1.*log(l1) - l1.*u - gammaln(h1));
f2 = exp( (h2-1).*log(u) +h2.*log(l2) - l2.*u - gammaln(h2))/p(5);
hrf = f1 - f2;
hrf = hrf((0:(p(7)/TR))*fMRI_T + 1); hrf = hrf'/sum(hrf);
hrf_plot_norm = hrf; % OK, this should be equivalent to spm defaults now.

%hrf_plot_norm_spm = spm_hrf(2); % and now go on from there
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
        task_events_file = ''; %#ok<NASGU> 
    return 
    end
end


% Now do things!
cd(output_dir)

% remove old ic_auto_selection folder if it exists
if isfolder('./ic_auto_selection')
    rmdir './ic_auto_selection' 's'
end

% Now do set up to make future coding easier and faster
% Describe the first (Factor) and second level layers
Factors = struct;
Layer_Factors = {'WholeBrain', 'InOutBrain', 'Networks', 'Regions', 'Power', 'Smoothness', 'Corr', 'Conditions'};
Layer_WholeBrain = {'fullvolume_smoothed_nothresh', 'fullvolume_nothresh'};
Layer_InOutBrain = {'Inbrain', 'Outbrain'};
Layer_Networks = {'MedialVisual', 'SensoryMotor', 'DorsalAttention', 'VentralAttention', 'FrontoParietal', 'DefaultModeNetwork', 'Subcortical'};
Layer_Regions = {'Signal', 'GM', 'WM', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'}; % all smaller regions
Layer_Power = {'Lowfreq', 'BOLDfreq', 'Highfreq'};
Layer_Smoothness = {'Smoothing_Retention'};
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

% (1) Check Spatial Map Overlap with ROIs
% (1a) Approximate overlap with mean and numvoxels

% We can just make tables for everything too
WholeBrain_labels = {'fullvolume_smoothed_nothresh', 'fullvolume_nothresh'};
ROI_labels = {'Signal', 'GM', 'WM', 'Edge', 'Subepe', 'CSF', 'Suscept', 'OutbrainOnly'};
InOut_labels = {'Inbrain', 'Outbrain'};
Freq_labels = {'Lowfreq', 'BOLDfreq', 'Highfreq'};
Networks_labels = {'MedialVisual', 'SensoryMotor', 'DorsalAttention', 'VentralAttention', 'FrontoParietal', 'DefaultModeNetwork', 'Subcortical'};

% Calculate ICsums for relevant things
ICsum_labels = [WholeBrain_labels, InOut_labels, Networks_labels, ROI_labels];
ICmean_table = table();
ICnumvoxels_table = table();
ICsum_table = table();
for i1 = 1:length(ICsum_labels)
    ICmean_table{:,ICsum_labels{i1}} = readmatrix(['ROIcalcs/', ICsum_labels{i1}, '_ICmean.txt']);
    ICnumvoxels_table{:,ICsum_labels{i1}} = readmatrix(['ROIcalcs/', ICsum_labels{i1}, '_ICnumvoxels.txt']);
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
% readmatrix some files first)
% Grab number of ICs and explained variance:
ICs = (1:length(ROI_genprops_table{:,'Signal'}))';
ICs_table = array2table(ICs, 'VariableNames', {'ICs'});
IC_exp_var = readmatrix('ROIcalcs/IC_exp_variance.txt');
IC_exp_var = IC_exp_var ./ 100; % put it in decimal places
IC_exp_var_table = array2table(IC_exp_var, 'VariableNames', {'Variance_Explained'});
% Add to Tables Struct
Tables.IC_table = ICs_table;
Tables.IC_exp_var_table = IC_exp_var_table;

% And also calculate smoothing values
Smoothness_table = table();
Smoothness_table{:, 'Smoothing_Retention'} = ICsum_table{:,'fullvolume_smoothed_nothresh'} ./ ICsum_table{:, 'fullvolume_nothresh'}; % this is the one we use
Tables.Smoothness_table = Smoothness_table;

% (2) Check Power Frequency Analysis in low, BOLD (0.008-0.15), and high
% frequencies
% (2a) Grab time series and calculate power and proportions in each
% frequency range of interest
N = T/TR; F = 1/TR;
df = 1/T; N_freq = N/2 + 1; %#ok<NASGU> 
f = F*(0:floor(N/2))/N; % to chart what frequencies we are at
frequencies = f; % easier naming later
lower_phys_cutoff = round(0.008 / df) + 1; % start at freq 0 at position 1
higher_phys_cutoff = round(0.15 / df) + 1;
cush = 1; %#ok<NASGU> % to avoid overlap, give an index of padding

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
P1_single_hrf_bp = (P1_single_hrf_bp ./ sum(P1_single_hrf_bp)); % renormalize sum based on bandpass - so it can be a fair comparison across other HRFs
Data.HRF_general.hrf_power_bp = P1_single_hrf_bp;
% Can plot(HRF_general.f,HRF_general.hrf_power_bp) to see the power if you want to confirm

% OK, now calculate power spectrums for each IC:
% set size of arrays explicitely to save computation time
ts = zeros(numvolumes, length(ICs));
Ps = zeros(floor(numvolumes/2)+1, length(ICs)); % Powerspectrums
f_all_power = zeros(length(ICs), 1);
BOLDfreqIC = zeros(length(ICs), 1);
LowfreqIC = zeros(length(ICs), 1);
HighfreqIC = zeros(length(ICs), 1);
for i=1:length(ICs)
    % grab time spectrum
    file=strcat('melodic/report/t', num2str(i), '.txt');
    ts(:,i)=readmatrix(file);

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
Data.N = N; Data.T = T; Data.tr = TR; Data.frequencies = f;
Data.lowfreq_cutoff = 0.008; Data.highfreq_cutoff = 0.15;

% Calculate overlap similarity and record in table, this may be added to
% later, consider not normalizing by range here for better outlier
% comparison later
HRF_table = table();
HRF_table{:,'general_power_overlap'} = sum(Data.HRF_general.hrf_power_bp .* Data.Powers)'; % a higher sum for an IC indicates more power spectrum overlap with expected BOLD response
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
Data.Confounds.Sixmotionparam = table2array(allconfounds(:,{'trans_x', 'trans_y', 'trans_z', 'rot_x', 'rot_y', 'rot_z'}));
Data.Confounds.DVARS = table2array(allconfounds(:,{'dvars'}));
Data.Confounds.FD = table2array(allconfounds(:,{'framewise_displacement'}));
Results.Confounds = Data.Confounds; %#ok<STRNU> 

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
condition_corrs_table = table();
if isfile(task_events_file)
    condition_names = unique(task_events_table{:,3})';
    block = zeros(numvolumes, length(condition_names));
    
    for j = 1:length(condition_names)
        curr_condition_name = condition_names{j};
        curr_condition_rows = find(strcmp(task_events_table{:,3}, curr_condition_name));
        for k = 1:length(curr_condition_rows)
            % ceiling is used incase it is some weird task that does not
            % fit TR lengths, round up to model a delayed HRF
            % start index add 1 because if start at timepiont 0, index is 1
            % in matlab
            start_ind = ceil(task_events_table{curr_condition_rows(k),1}./TR)+1;
            % don't add one for stop index, because in relation to start
            % index
            stop_ind = start_ind + ceil(task_events_table{curr_condition_rows(k),2}./TR);
            block(start_ind:stop_ind, j) = 1;
        end
    end

    events_table = array2table(block, 'VariableNames', condition_names);

    % The above should give us blocks for all conditions, now convolve with hrf
    % to get a guestimated hrf response for each condition:
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
    ICs_reps = repmat(ICs, 1, size(condition_corrs,2)); %#ok<NASGU> 
    
    % both calculate relative correlation strength for hrfs AND calculate
    % power spectra
    % initialize things first, then go into for loop
    N = T/TR; F = 1/TR; df = 1/T; N_freq = N/2 + 1; %#ok<NASGU> 
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
        P1_all_hrfs_bp(:,m) = (P1_all_hrfs_bp(:,m) ./ sum(P1_all_hrfs_bp(:,m))); % renormalize sum based on bandpass for fair future hrf comparisons
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
    [~,I] = max(condition_corrs, [], 2);
    condition_corr_most_correlated_index = I;
    max_conditioncorr_powerspectra = P1_all_hrfs_bp(:,I)'; % each row is an IC, each column is a positive powerspectrum

    % OK, NOW we can calculate relative poweroverlap, consider not
    % normalizing "range" for better comparison across subjects later
    % idea is helps us find outliers in subjects later, where power overlap
    % is simply more poor
    power_overlap_all_hrfs = sum(max_conditioncorr_powerspectra .* Ps', 2); % a higher sum for an IC indicates more power spectrum overlap with expected BOLD response
    hrfstask_powerspectrum = P1_all_hrfs_norm;
    hrfstask_fullpowerspectrum = P2_all_hrfs;
    HRF_table{:, 'besttask_power_overlap'} = power_overlap_all_hrfs;

    
    save('Task_Condition_Info.mat', 'ts', 'numvolumes','hrf_conditions_table', 'condition_corr_most_correlated_index', ...
        'condition_corrs', 'condition_names', 'hrfstask_fullpowerspectrum', ...
        'hrfstask_powerspectrum', 'max_conditioncorr_powerspectra', 'power_overlap_all_hrfs', 'frequencies', 'hrf_conditions')
end

% save standard 9Param and 8Param for comparison or use later, and write to mat file for fsl glm in future:
writematrix(Data.Confounds.Nineparam, '9p_regressors.txt', 'Delimiter', ' ')
writematrix(Data.Confounds.Eightparam, '8p_regressors.txt', 'Delimiter', ' ')
writematrix(Data.Confounds.Sixmotionparam, '6p_motion_regressors.txt', 'Delimiter', ' ')
% include versions with an intercept column, in case you regress without
% demeaning 
intercept = ones(size(Data.Confounds.Eightparam,1),1);
writematrix([intercept, Data.Confounds.Nineparam], '9p_regressors_intercept.txt', 'Delimiter', ' ')
writematrix([intercept, Data.Confounds.Eightparam], '8p_regressors_intercept.txt', 'Delimiter', ' ')
writematrix([intercept, Data.Confounds.Sixmotionparam], '6p_motion_regressors_intercept.txt', 'Delimiter', ' ')
[~, ~] = call_fsl('Text2Vest 8p_regressors_intercept.txt 8p_regressors_intercept.mat');
[~, ~] = call_fsl('Text2Vest 9p_regressors_intercept.txt 9p_regressors_intercept.mat');


% Calculate the highest normalized HRF matching (whether to a general HRF
% frequency range, or HRFs in response to tasks)
% get a nromalized range of HRF_table first
HRF_table_norm = normalize(Tables.HRF_table, "range");
HRF_table{:, 'best_power_overlap_norm'} = max(HRF_table_norm{:,:}, [], 2); % if there is no task, it will just be the same as general but normalized
HRF_table_norm{:, 'best_power_overlap_norm'} = HRF_table.best_power_overlap_norm;
Tables.HRF_table = HRF_table;
Tables.HRF_table_norm = HRF_table_norm;


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
    max_feats = find(val_C == max(val_C)); % make separate variables to deal with the unlikely scenario 3 groups becomes only 2. Likely terrible data if that happens anyways.
    min_feats = find(val_C == min(val_C));
    classification_table{:, ['High_', feature_labels{i1}]} = val_idx == max_feats(end);
    classification_table{:, ['Low_', feature_labels{i1}]} = val_idx == min_feats(1);
end

% Start parsing out what is important to look at from classification
best_classification_labels = {'High_Signal', 'High_Smoothing_Retention', 'High_best_power_overlap_norm'};
% there is no best task power if there are no condition cors
if ~isempty(condition_corrs)
    othergood_classification_labels = {'High_GM', 'High_general_power_overlap', 'High_besttask_power_overlap', 'High_WM'};
else
    othergood_classification_labels = {'High_GM', 'High_general_power_overlap', 'High_WM'};
end

bad_classification_labels = {'Low_Signal', 'Low_GM', 'High_Edge', 'High_Subepe', 'High_CSF', 'High_Suscept', ...
    'High_OutbrainOnly', 'High_Outbrain', 'Low_BOLDfreq', 'High_Highfreq', ...
    'Low_best_power_overlap_norm', 'Low_Smoothing_Retention', 'High_DVARS_Corr', 'High_FD_Corr'}; 
networks_classification_labels = {'High_MedialVisual', 'High_SensoryMotor', ...
    'High_DorsalAttention', 'High_VentralAttention', 'High_FrontoParietal', ...
    'High_DefaultModeNetwork', 'High_Subcortical'};

condition_corrs_classification_labels = cell(1, length(condition_names));
for i1 = 1:length(condition_names)
    condition_corrs_classification_labels{i1} = ['High_', condition_names{i1}];
end

best_classification_table = classification_table(:, best_classification_labels);
othergood_classification_table = classification_table(:, othergood_classification_labels);
good_classification_table = [best_classification_table, othergood_classification_table];
bad_classification_table = classification_table(:, bad_classification_labels);
networks_classification_table = classification_table(:, networks_classification_labels);
% There are no condition corrs classification if there are no condition
% corrs
if ~isempty(condition_corrs)
    condition_corrs_classification_table = classification_table(:, condition_corrs_classification_labels);
else
    condition_corrs_classification_table = table(); % no condition corrs, make table empty
end


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

% check if the only "bad thing" is high_Lowfreq - this alone is not
% worrisome
%not_so_bad_indices = (sum(table2array(bad_classification_table), 2) == 1) & bad_classification_table.High_LowFreq

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
[~, signal_idx] = sort(normalize(feature_relative_table{:,'Smoothing_Retention'}, "range") .* ...
    normalize(feature_relative_table{:,'Signal'}, "range").^2 .* ...
    normalize(feature_relative_table{:,'best_power_overlap_norm'}, "range"), 'descend');

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
    & (classification_final_table{:, 'High_best_power_overlap_norm'} == 1);
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
            signal_decision(g) = 0; %#ok<AGROW> 
        else
            stop_num = stop_num + 1;
            % mark as signal
            signal_decision(g) = 1; %#ok<AGROW> 
        end
    elseif ismember(signal_idx(g), HasSig_ICs) == 1
        % Needs to have at least one signal-like marking
        stop_num = stop_num + 1;
        % mark as signal
        signal_decision(g) = 1; %#ok<AGROW> 
    else
        % if it has no signal-like markings, mark as noise
        stop_num = stop_num - 1;
        % mark as signal
        signal_decision(g) = 0; %#ok<AGROW> 
    end
       % if stop num went above tolerance, set it back to tolerance
       if stop_num > tolerance
           stop_num = tolerance;
       end
    g = g + 1; % g will be how far along the signal ICs you go
end
g = g - 1;

potential_signal_ICs = signal_idx(1:g);
signal_indices = zeros(1,length(ICs));
signal_indices(potential_signal_ICs) = signal_decision; % Labels based on signal decision earlier, and anything past that as noise
signal_indices = logical(signal_indices');


% Give a fail safe incase it is so bad that no ICs are labeled as signal
if sum(signal_indices) == 0
    fprintf('WARNING: No ICs were labeled as signal. Will label first IC in best order as signal to allow the rest of the code to run.\n')
    signal_indices(:,1) = true;
end
noise_indices = ~signal_indices;
signal_ICs = ICs(signal_indices);
noise_ICs = ICs(noise_indices);

signal_checker_table = array2table(signal_idx, 'VariableNames', {'PotentialICs'});
signal_checker_table.SignalLabel = signal_indices(signal_idx);
signal_checker_table.HighSignalLabel = HighSig_indices(signal_idx);
signal_checker_table.HighNoiseLabel = Highnoise_indices(signal_idx);
Highconditioncorr_indices = sum(table2array(condition_corrs_classification_table),2) > 0;
if ~isempty(condition_names)
    signal_checker_table.TaskCorrelated = double(sum(Highconditioncorr_indices(signal_idx,:),2) > 0);
end


% add on classification labels to signal_checker_table 
% One column for good, one for bad, one for networks, one for tasks (if it
% exists)
% good first:
good_classification_potentialsignal_labels = cell(length(signal_idx), 1);
for n = 1:length(signal_idx)
    good_class_labels = logical(good_classification_table{signal_idx(n),:});
    % put 'None' if none exist, helps visualize on excel
    if sum(good_class_labels) == 0
        curr_class = {'None'};
    else
        %curr_class = {good_classification_table.Properties.VariableNames{good_class_labels}};
        curr_class = good_classification_table.Properties.VariableNames(good_class_labels);
    end
    %good_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
    good_classification_potentialsignal_labels{n} = join(curr_class(:), '; ');
end
signal_checker_table.Good_Tags = good_classification_potentialsignal_labels;

% now bad:
bad_classification_potentialsignal_labels = cell(length(signal_idx), 1);
for n = 1:length(signal_idx)
    bad_class_labels = logical(bad_classification_table{signal_idx(n),:});
    % put 'None' if none exist, helps visualize on excel
    if sum(bad_class_labels) == 0
        curr_class = {'None'};
    else
        %curr_class = {bad_classification_table.Properties.VariableNames{bad_class_labels}};
        curr_class = bad_classification_table.Properties.VariableNames(bad_class_labels);
    end
    %bad_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
    bad_classification_potentialsignal_labels{n} = join(curr_class(:), '; ');
end
signal_checker_table.Bad_Tags = bad_classification_potentialsignal_labels;

% and networks
networks_classification_potentialsignal_labels = cell(length(signal_idx), 1);
for n = 1:length(signal_idx)
    networks_class_labels = logical(networks_classification_table{signal_idx(n),:});
    % put 'None' if none exist, helps visualize on excel
    if sum(networks_class_labels) == 0
        curr_class = {'None'};
    else
        %curr_class = {networks_classification_table.Properties.VariableNames{networks_class_labels}};
        curr_class = networks_classification_table.Properties.VariableNames(networks_class_labels);
    end
    %networks_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
    networks_classification_potentialsignal_labels{n} = join(curr_class(:), '; ');
end
signal_checker_table.Networks_Tags = networks_classification_potentialsignal_labels;

% finally condition corrs (if they exist)
if ~isempty(condition_names)
    condition_corrs_classification_potentialsignal_labels = cell(length(signal_idx), 1);
    for n = 1:length(signal_idx)
        condition_corrs_class_labels = logical(condition_corrs_classification_table{signal_idx(n),:});
        % put 'None' if none exist, helps visualize on excel
        if sum(condition_corrs_class_labels) == 0
            curr_class = {'None'};
        else
            %curr_class = {condition_corrs_classification_table.Properties.VariableNames{condition_corrs_class_labels}};
            curr_class = condition_corrs_classification_table.Properties.VariableNames(condition_corrs_class_labels);
        end
        %condition_corrs_classification_potentialsignal_labels{n} = join({curr_class{:}}, '; ');
        condition_corrs_classification_potentialsignal_labels{n} = join(curr_class(:), '; ');
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

% estimate effective degrees of freedom (somewhat) from CICADA alone via proportions:
Results.percent_variance_kept = sum(IC_exp_var(signal_ICs)); % probably want > 0.1 (10%)
Results.num_ICs_kept = length(signal_ICs);
Results.num_ICs_total = length(ICs);
Results.IC_exp_var = IC_exp_var;
Results.ICs = ICs;

% Create images to review signal and noise components from
% probabilities ICs
all_prob = niftiread('./melodic/ICprobabilities.nii.gz');
gm_prob = niftiread('./region_masks/GM_prob.nii.gz');
notgm_prob = niftiread('./region_masks/NotGM_prob.nii.gz');
gm_bin = gm_prob > 0.75;
notgm_bin = notgm_prob > 0.75;
noise_prob = all_prob(:,:,:,noise_ICs);
noise_prob_info = niftiinfo('./melodic/ICprobabilities.nii.gz');
signal_prob = all_prob(:,:,:,signal_ICs);

% grab noise probabilities, update header, and write, same for
% signal
noise_prob_info.ImageSize = size(noise_prob);
noise_prob_info.DisplayIntensityRange = [0.9 1];
noise_prob_info.Datatype = 'single';

signal_prob_info = noise_prob_info;
signal_prob_info.ImageSize = size(signal_prob);

potential_signal_prob = all_prob(:,:,:,signal_idx);
potential_signal_prob_info = noise_prob_info;
potential_signal_prob_info.ImageSize = size(potential_signal_prob);

niftiwrite(potential_signal_prob, 'PotentialSignalICs', potential_signal_prob_info, 'Compressed', true) % just the ICs that were considered
niftiwrite(noise_prob, 'NoiseICs', noise_prob_info, 'Compressed', true)
niftiwrite(signal_prob, 'SignalICs', signal_prob_info, 'Compressed', true)

% Now compress noise and signal, with only max probabilities per voxel
[~, ~] = max(all_prob, [], 4); 
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
writetable(compare_cleaning, 'compare_auto_cleaning.csv', 'WriteRowNames', true)

% and then make noise components if you want to do aggressive denoising like in CONN (regression)
mixing_matrix = readmatrix('./melodic/melodic_mix');
auto_noise_dist_covariates = mixing_matrix(:, Results.noise_ICs');
writematrix(auto_noise_dist_covariates, 'auto_noise_dist_covariates.txt', 'Delimiter', ' ')
writematrix(auto_noise_dist_covariates, 'auto_noise_dist_covariates.csv')

movefile auto* ic_auto_selection
movefile *SignalIC* ic_auto_selection
movefile *NoiseIC* ic_auto_selection
movefile *IC_auto_checker.* ic_auto_selection
movefile *_auto_* ic_auto_selection
movefile *regressors* regressors_timeseries
movefile *timeseries.* regressors_timeseries
% move confounds back out though so it is easier to find
cd regressors_timeseries
movefile confounds_timeseries* ../
cd ../

% save a matrix of relevant variables if you want to examine later!
save('DecisionVariables_Auto.mat', 'output_dir', 'subject_id', 'session_id', 'task_id', 'ICs', 'Data', 'Tables', 'Results')

movefile *Auto.* ic_auto_selection

% OK, now run the denoising: Start with CICADA
task_dir = pwd;
prefix = [subject_id, '_', session_id, '_task-', task_id];
CICADA_tag = 'CICADA_auto_nonagg';
suffix = 'bold.nii.gz';
fsl_regfilt_command = ['fsl_regfilt -i ', task_dir, '/funcfile.nii.gz -f ', ...
    '"', '$(cat ', task_dir, '/ic_auto_selection/auto_noise_dist_ICs.csv)', '"', ...
    ' -d ', task_dir, '/melodic/melodic_mix -m ', ...
    task_dir, '/funcmask.nii.gz -o ', prefix, '_', CICADA_tag, '_', suffix];
fprintf(['Running: ', fsl_regfilt_command, '\n'])
[~, ~] = call_fsl(fsl_regfilt_command);

% grab Tmean to add back in
tmean_command = ['fslmaths ', task_dir, '/funcfile.nii.gz -Tmean ', task_dir, '/tmean_funcfile.nii.gz'];
[~, ~] = call_fsl(tmean_command);

% Now do 8 parameter & 9 parameter
eightparam_regress_command = ['fsl_glm -i ', task_dir, '/funcfile.nii.gz -d ', task_dir, ...
    '/regressors_timeseries/8p_regressors_intercept.mat -m ', task_dir, '/funcmask.nii.gz ', ...
    '--out_res=', prefix, '_8p_', suffix];
fprintf(['Running: ', eightparam_regress_command, '\n'])
[~, ~] = call_fsl(eightparam_regress_command);
nineparam_regress_command = ['fsl_glm -i ', task_dir, '/funcfile.nii.gz -d ', task_dir, ...
    '/regressors_timeseries/9p_regressors_intercept.mat -m ', task_dir, '/funcmask.nii.gz ', ...
    '--out_res=', prefix, '_9p_', suffix];
fprintf(['Running: ', nineparam_regress_command, '\n'])
[~, ~] = call_fsl(nineparam_regress_command);

% add Tmean back onto 8 & 9 parameter 
tmean_add_8p_command = ['fslmaths ', prefix, '_8p_', suffix, ' -add ', task_dir, '/tmean_funcfile.nii.gz ', prefix, '_8p_', suffix];
tmean_add_9p_command = ['fslmaths ', prefix, '_9p_', suffix, ' -add ', task_dir, '/tmean_funcfile.nii.gz ', prefix, '_9p_', suffix];
[~, ~] = call_fsl(tmean_add_8p_command);
[~, ~] = call_fsl(tmean_add_9p_command);

% create cleaned directory
cleaned_dir = [output_dir, '/cleaned'];
if ~isfolder(cleaned_dir)
    mkdir(cleaned_dir)
end
delete([cleaned_dir, '/*']) % clear out old files in here, save space and just in case

% for good measure, put in a copy of the original with similar naming in cleaned dir - Good for comparisons in QC script
copyfile('funcfile.nii.gz', [cleaned_dir, '/', prefix, '_orig_', suffix])
movefile('*CICADA*.nii.gz', cleaned_dir)
movefile('*8p*.nii.gz', cleaned_dir)
movefile('*9p*.nii.gz', cleaned_dir)

% get cleaned_filename
cleaned_file_info = dir([cleaned_dir, '/*CICADA*auto*.nii.gz']);
cleaned_file = [cleaned_file_info.folder, '/', cleaned_file_info.name];

end