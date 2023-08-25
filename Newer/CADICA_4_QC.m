% CADICA_4_QC
% CADICA_3_Clean calculated parcellation text files from Gordon at last of
% denoised files. This script just reads them and runs correlations to plot
% QC

clearvars
home='/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA_Updated';
subjects = {'102'};
sessions = {'01'};

cd(home)

for k = 1:length(subjects)
    cd(home)
    currsubjfol = ['sub-', subjects{k}];
    cd(currsubjfol)

    for l = 1:length(sessions)
        currsessfol = ['ses-', sessions{l}];
        cd(currsessfol)

       % confounds:
        confound_place = './confounds_timeseries.csv';
        allconfounds = readtable(confound_place);
        confounds_dvars = table2array(allconfounds(:,{'dvars'}));
        
        % Get GM prob mask and create GM mask and not GM mask for denoised
        % and original files
        denoised_file = './cleaned/s_bp_sub-102_ses-01_task-rest_space-MNI152NLin2009cAsym_desc-CADICA_agg_bold.nii.gz';
        GM_prob_file = './anatmasks/GMprob_resam.nii.gz';
        funcmask = './funcmask.nii.gz';
        [denoised_GM, denoised_NotGM, denoised_GS] = getData(denoised_file, funcmask, GM_prob_file);

        orig_file = '../funcfile.nii.gz';
        [orig_GM, orig_NotGM, orig_GS] = getData(orig_file, funcmask, GM_prob_file);
        
        % Now compute relevant correlations
        [denoised_GM_NotGM_corr, denoised_GM_GM_corr, denoised_GM_dvars_corr] = createHistData(denoised_GM, denoised_NotGM, dvars, 500);
        [orig_GM_NotGM_corr, orig_GM_GM_corr, orig_GM_dvars_corr] = createHistData(orig_GM, orig_NotGM, dvars, 500);
        

        % GM_NotGM corr: Want to see shift toward center in denoised
        % (Remove connections between GM and not GM - noise)
        figure
        hold on
        title('Original vs Denoised GM-NotGM Corrs')
        histogram(orig_GM_NotGM_corr, 'Normalization', 'pdf')
        histogram(denoised_GM_NotGM_corr, 'Normalization', 'pdf')
        legend('Original', 'Denoised')
        hold off

        % GM_GM corr: Do not want to see a complete shift toward center in
        % denoised (Retain true GM signal connections)
        figure
        hold on
        title('Original vs Denoised GM-GM Corrs')
        histogram(orig_GM_GM_corr, 'Normalization', 'pdf')
        histogram(denoised_GM_GM_corr, 'Normalization', 'pdf')
        legend('Original', 'Denoised')
        hold off
        
        % GM_dvars corr: Want to see a shift toward center for denoised
        % (removal of noise/sudden GS signal shift in GM)
        figure
        hold on
        title('Original vs Denoised GM-DVARS Corrs')
        histogram(orig_GM_dvars_corr, 'Normalization', 'pdf')
        histogram(denoised_GM_dvars_corr, 'Normalization', 'pdf')
        legend('Original', 'Denoised')
        hold off

        % Finally, Plot GS before and after Denoising (spikes should be
        % greatly diminished
        figure
        hold on
        title('Original vs Denoised Global Signal')
        plot(denoised_GS)
        plot(orig_GS)
        legend('Original', 'Denoised')
        hold off
        
    end
end

function [GM, NotGM, Global_Signal] = getData(funcfile, funcmask, GMprobfile)
     % read in the niftis
     GM_prob_data = niftiread(GMprobfile);
     func_data = niftiread(funcfile);
     funcmask_data = niftiread(funcmask);

    GM_mask = logical(GM_prob_data > 0.67 & funcmask_data == 1);
    NotGM_mask = logical(GM_prob_data < 0.33 & funcmask_data == 1);
    
    % Select the timeseries for GM and not GM and global signal
    global_signal = zeros(sum(funcmask_data(:)), size(func_data, 4));
    GM = zeros(sum(GM_mask(:)), size(func_data, 4));
    NotGM = zeros(sum(NotGM_mask(:)), size(func_data, 4));
    for j = 1:size(func_data, 4)
        curr_func = func_data(:,:,:,j);
        GM(:,j) = curr_func(GM_mask);
        NotGM(:,j) = curr_func(NotGM_mask);
        global_signal(:,j) = curr_func(funcmask_data == 1);
    end
    
    GM = GM(GM(:,1) ~= 0, :);
    NotGM = NotGM(NotGM(:, 1) ~= 0, :);
    Global_Signal = global_signal(global_signal(:, 1) ~= 0, :);

end

function [GM_NotGM_corr, GM_GM_corr, GM_dvars_corr] = createHistData(GM, NotGM, dvars, perms)
        % randperm selections
        GM_randperm = randperm(size(GM, 1));
        NotGM_randperm = randperm(size(NotGM, 1));

        % Should center close to 0 (not globally correlated) if it's
        % denoised (e.g., denoised lowered noise)
        GM_NotGM_corr = corr(GM(GM_randperm(1:perms), :)', NotGM(NotGM_randperm(1:perms), :)');
        GM_NotGM_corr = GM_NotGM_corr(tril(GM_NotGM_corr, -1) ~= 0);

        % Should center above 0 (e.g., denoise did not remove true signal)
        GM_corr = corr(GM(GM_randperm(1:perms), :)', GM(GM_randperm(1:perms), :)');
        GM_GM_corr = GM_corr(tril(GM_corr, -1) ~= 0);

        % corr dvars vs GM diff abs: Should be reduced if denoised
        GM_dvars_corr = corr(dvars(2:end), abs(diff(GM(GM_randperm(1:perms), :)')));
        GM_dvars_corr = GM_dvars_corr(tril(GM_corr, -1) ~= 0);
end