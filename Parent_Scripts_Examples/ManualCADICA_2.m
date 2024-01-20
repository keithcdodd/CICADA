% Running CADICA Denoising Part 2 of 2 (runs through manual cleaning,
% following auto cleaning)
% 
% NOTE: NEED TO MAKE IC_manual_checker.csv BEFORE RUNNING THIS SCRIPT:
% To make IC_manual_checker.csv, open IC_auto_checker.csv in
% ic_auto_selection folder (assuming you have run the auto version already). 
% Use melodic report, IC_probabilities, and any
% other information you might want to make decisions to edit the Potential
% Signal Labels column in IC_auto_checker.csv (1 means signal, 0 is noise) 
% [we suggest using Griffanti et al. 2017 Hand classification of fMRI ICA 
% noise components paper as a reference], and then save as 
% IC_manual_checker.csv in that same location (ic_auto_selection folder). 
% THEN you can run this script.

clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%
subj_ids = {'122'};
%subj_ids = {'144', '181', '184', '185', '186'};
sess_ids = {'01'};
%sess_ids = {'01', '02', '03', '04'}; %, '02', '03', '04'};
task_names = {'rest'}; % include runs if they exist, e.g., 'foodpics_run-01'
tr = 2;
cadica_dir = '/data/images/awesome/data/bids_data/derivatives/cadica_2mm_stc_kd'; % 2mm resolution output space, slice time corrected, Keith Dodd
% Decide which scripts you want to run.
run_CADICA_2B_ManualLabeling = 1; % 0 if you don't want to run it, e.g., it was already run
run_CADICA_3_ManualClean = 1;
run_CADICA_4_QC = 1;
lowfreq_cutoff = 0.008; % Hz, 0.008 default, relevant for 3rd and 4th function
highfreq_cutoff = 0.15; % Hz, 0.15 default, relevant for 3rd and fourth function
blur = 6; % 6mm gaussian, 6mm default, relevant for 3rd and fourth function
ic_filter_type = 'nonagg'; % nonagg or agg for comparisons, only relevant for 4th function
filterchoice = 'nf'; % nf is no filter, hp is highpass, bp is bandpass. All are calculated, this is just for comparisons for 4_QC, relevant for 4th function
smoothchoice = 'ns'; % ns (no smooth) or s (smooth) for comparison for QC, relevant for fourth function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% For FSL to run all happy, we need to set the LD_LIBRARY_PATH Differently.
% But just for this script.

setenv('LD_LIBRARY_PATH', '/usr/lib64/lib64/libstdc++.so.6:/usr/local/MATLAB/R2022b/bin/matlab')

CADICA_3_path = '/home/doddke/code/CADICA/CADICA_3_Clean.sh';
addpath('/home/doddke/code/CADICA')


for j = 1:length(subj_ids)

    fprintf('\n')
    fprintf(['Running for subject ', subj_ids{j}, '\n'])
    % check for folder
    curr_subj_fol = [cadica_dir, '/sub-', subj_ids{j}];
    if ~isfolder(curr_subj_fol)
        fprintf(['No subject folder found for ', subj_ids{j}, '. Moving on to next one...\n'])
        continue
    end

    for k = 1:length(sess_ids)
        
        fprintf(['Running for session ', sess_ids{k}, '\n'])
        % check for folder
        curr_sess_fol = [curr_subj_fol, '/ses-', sess_ids{k}];
        if ~isfolder(curr_sess_fol)
            fprintf(['No sess folder found for ', sess_ids{k}, '. Moving on to next one...\n'])
            continue
        end

        for l = 1:length(task_names)
            fprintf(['Running for task ', task_names{l}, '\n'])

            curr_dir = [curr_sess_fol, '/', task_names{l}];

            % First check if ic_manual_checker.csv exists, because if not,
            % skip and go to next
            IC_manual_checker = [curr_dir, '/ic_auto_selection/IC_manual_checker.csv'];
            if ~isfile(IC_manual_checker)
                IC_manual_checker = [curr_dir, '/ic_manual_selection/IC_manual_checker.csv'];
                if ~isfile(IC_manual_checker)
                    fprintf('Cannot find IC_manual_checker. Continuing to next one...\n')
                    continue
                end
            end

           
            % For Manual, you do not rerun the first script, go right to
            % the modified manual second script
            
            
            % Now run the Manual labeling
            if run_CADICA_2B_ManualLabeling == 1
                fprintf('Running: CADICA_2B_ManualLabeling\n')
                CADICA_2B_ManualLabeling(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l})
                fprintf('Done with CADICA_2B_ManualLabeling\n')
            end
            
            % Now run the Cleaning of the Manual labeling
            if run_CADICA_3_ManualClean == 1
               CADICA_3_command = [CADICA_3_path, ' -c ', cadica_dir, ' -p ', subj_ids{j}, ' -s ', sess_ids{k}, ' -t ', task_names{l}, ' -i ', 'manual', ' -l ', num2str(lowfreq_cutoff), ' -h ', num2str(highfreq_cutoff), ' -b ', num2str(blur), ' < /dev/null'];
               fprintf(['Running: ', CADICA_3_command, '\n'])
               [status, cmdout_CADICA_3] = system(CADICA_3_command, '-echo');
            end
            
            if run_CADICA_4_QC == 1
               fprintf('Running: CADICA_4_QC\n')
               % NOTE: [prefix, compare_file_tag, suffix] should be the
               % full string for the compare filename
               prefix = ['sub-', subj_ids{j}, '_ses-', sess_ids{k}, '_task-', task_names{l}, '_'];
               suffix ='_bold.nii.gz';

               % Compare to 8p first (6 movement, CSF, and WM)
               CADICA_4_QC(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, prefix, suffix, '8p', 'manual', ic_filter_type, filterchoice, smoothchoice, lowfreq_cutoff, highfreq_cutoff, tr);
               
               % also compare to 9p (with global signal)
               CADICA_4_QC(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, prefix, suffix, '9p', 'manual', ic_filter_type, filterchoice, smoothchoice, lowfreq_cutoff, highfreq_cutoff, tr);
               fprintf('Done with CADICA_4\n\n')

               % then compare to the auto!
               auto_tag = ['CADICA_auto_', ic_filter_type];
               CADICA_4_QC(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, prefix, suffix, auto_tag, 'manual', ic_filter_type, filterchoice, smoothchoice, lowfreq_cutoff, highfreq_cutoff, tr);
               
            end
            fprintf('Done with task ', task_names{l}, ' for session ', sess_ids{k}, ' for subject ', subj_ids{j}, '\n')
            fprintf('\n\n\n')
            
        end
    end
end


