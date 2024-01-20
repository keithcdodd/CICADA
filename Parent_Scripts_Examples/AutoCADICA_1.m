% Running CADICA Denoising Part 1 of 2 (runs through autolabeling)

% After this is done, go to Part 2 to see instructions for manual labeling,
% and then running the manual labeling edits, alongside cleaning, and QC
clearvars
%%%%%%%%%%%%%%%%%%%%%%%%%% SET UP %%%%%%%%%%%%%%%%%%%%%%%%%%
fmriprep_dir = '/data/images/awesome/data/bids_data/derivatives/fmriprep_2mm_stc_kd'; % 2mm resolution output space, slice time corrected, Keith Dodd
subj_ids = {'102', '103', '105', '106', '107', '108', '109', '112', '113', '114', '115', '116', '117', '118', '121', '122'};
%subj_ids = {'164', '165', '168', '169', '171', '172', '173', '174', '176', '178', '181', '184', '185', '186', '187'};
%subj_ids = {'122'};
sess_ids = {'01', '02', '03', '04'};
%sess_ids = {'01'};
task_names = {'foodpics_run-02'}; % include runs if they exist, e.g., 'foodpics_run-01', don't put both runs at once though!
%numvolumes = 300; % WILL LIKELY DIFFER DEPENDING ON TASK: how many timepoints for task (this is why you would likely not put both rest and task data together, as their timepoints would likely differ)
numvolumes = 172;
% task events file should be the one that you put into for fMRIprep (same
% structure)
task_events_file = '/data/images/awesome/data/bids_data/task-foodpics_run-01_events.tsv'; %'/data/images/awesome/data/bids_data/task-foodpics_run-01_events.tsv'; % put '' if you aren't using a task file. It is not necessary, just potentially helpful.
%task_events_file = '';
% for task_events_table_power: Make sure the names exactly match your task
% events file, the power values (0-1) allow for relative predictions of response.
% At very least do 0 for baseline, and 1 for the main task
% if there is a baseline condition, make sure the task_events_file has it
% labeled exactly as 'baseline' in column 3 (so we can remove it easily)
CADICA_script_home = '/home/doddke/code/CADICA'; % where the scripts are held
tr = 2; % repetition time
template_space = 'MNI152NLin6Asym'; % this needs to match the naming convention from fmriprep, a likely default is MNI152NLine6Asym
resolution = 2; % mm resolution of output space for the files you will use
anat_sess = '01'; % which session ID holds the best t1 anatomical to use as reference
cadica_dir = '/data/images/awesome/data/bids_data/derivatives/cadica_2mm_stc_kd'; % 2mm resolution output space, slice time corrected, Keith Dodd
redo_melodic = 0; % if 0, don't redo MELODIC if it has been done before, if 1, then redo it. Relevant for function 1.
num_ICs = 0; % 0 is default to do automatic estimation of ICs to calculate. You can give it a positive integer instead to tell MELODIC how many ICs to calculate.
tolerance = 5; % integer: higher numbers gives you more ICs you want to examine, lower number gives you less, but might miss some good ICs. 5 or 6 is pretty good. Relevant for function 2.
% Decide which scripts you want to run. Likely you want to run the first
% two scripts, but this will allow to just run one of them.
run_CADICA_1_ICs = 1; % 0 if you don't want to run it
run_CADICA_2_AutoLabeling = 1; % 0 if you don't want to run it
run_CADICA_3_AutoClean = 1;
run_CADICA_4_QC = 1;
lowfreq_cutoff = 0.008; % Hz, 0.008 default, only relevant for functions 3 and 4
highfreq_cutoff = 0.15; % Hz, 0.15 default, only relevant for functions 3 and 4
blur = 6; % 6mm gaussian, 6mm default, only relevant for functions 3 and 4
ic_filter_type = 'nonagg'; % nonagg or agg, nonagg default, only relevant for function 4
filterchoice = 'nf'; % nf is no filter, hp is highpass, bp is bandpass. All are calculated, this is just for comparisons for 4_QC, only relevant for function 4
smoothchoice = 'ns'; % ns (no smooth) or s (smooth) for comparison for QC, only relevant for function 4
% For FSL to run all happy, we need to set the LD_LIBRARY_PATH Differently.
% But just for this script (you will likely need to adjust this to your
% system, you can google FSL LD_LIBRARY_PATH Matlab to find more
% information)
setenv('LD_LIBRARY_PATH', '/usr/lib64/lib64/libstdc++.so.6:/usr/local/MATLAB/R2022b/bin/matlab')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CADICA_1_path = [CADICA_script_home, '/CADICA_1_MasksandICAs.sh'];
CADICA_3_path = [CADICA_script_home, '/CADICA_3_Clean.sh'];
addpath(CADICA_script_home)


for j = 1:length(subj_ids)
    fprintf('\n')
    fprintf(['Running for subject ', subj_ids{j}, '\n'])
    % check for folder
    curr_subj_fol = [fmriprep_dir, '/sub-', subj_ids{j}];
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
            fprintf('\n')

            % ready the command based on if you will redo MELODIC or not
            if redo_melodic == 0
                if num_ICs < 1
                    % if < 1, do automatic estimation for number of ICs
                    CADICA_1_command = [CADICA_1_path, ' -f ', fmriprep_dir, ' -p ',  subj_ids{j}, ' -s ', sess_ids{k}, ' -t ', task_names{l}, ' -r ', num2str(tr), ' -l ', ...
                    template_space, ' -d ', num2str(resolution), ' -a ', anat_sess, ' -o ', cadica_dir, ' < /dev/null'];
                else
                    % otherwise estimate the number of ICs given
                    CADICA_1_command = [CADICA_1_path, ' -f ', fmriprep_dir, ' -p ',  subj_ids{j}, ' -s ', sess_ids{k}, ' -t ', task_names{l}, ' -r ', num2str(tr), ' -l ', ...
                    template_space, ' -d ', num2str(resolution), ' -a ', anat_sess, ' -o ', cadica_dir, '-n ', num_ICs, '< /dev/null'];
                end
                
                
            elseif redo_melodic == 1
                if num_ICs < 1
                    % if < 1, do automatic estimation for number of ICs
                    CADICA_1_command = [CADICA_1_path, ' -f ', fmriprep_dir, ' -p ',  subj_ids{j}, ' -s ', sess_ids{k}, ' -t ', task_names{l}, ' -r ', num2str(tr), ' -l ', ...
                    template_space, ' -d ' num2str(resolution), ' -a ', anat_sess, ' -o ', cadica_dir, ' -m ', '< /dev/null'];
                else
                    % otherwise estimate the number of ICs given
                    CADICA_1_command = [CADICA_1_path, ' -f ', fmriprep_dir, ' -p ',  subj_ids{j}, ' -s ', sess_ids{k}, ' -t ', task_names{l}, ' -r ', num2str(tr), ' -l ', ...
                    template_space, ' -d ', num2str(resolution), ' -a ', anat_sess, ' -o ', cadica_dir, '-n ', num_ICs, ' -m ', '< /dev/null'];
                end
                
            end
            % Now you can run the first script!
            if run_CADICA_1_ICs == 1
                fprintf(['Running: ', CADICA_1_command, '\n'])
                [status, cmdout_CADICA_1] = system(CADICA_1_command, '-echo');
                fprintf('Done with CADICA_1\n\n')
            end
            
            
            % Now run the autolabeling
            if run_CADICA_2_AutoLabeling == 1
                if isfile(task_events_file)
                    fprintf('Running: CADICA_2_AutoLabeling with Task Events File\n')
                    CADICA_2_AutoLabeling(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, numvolumes, tr, tolerance, task_events_file)
                    fprintf('Done with CADICA_2\n\n')
                else
                    fprintf('Running: CADICA_2_AutoLabeling without Task Events File\n')
                    CADICA_2_AutoLabeling(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, numvolumes, tr, tolerance)
                    fprintf('Done with CADICA_2\n\n')
                end
            end
            
            % Now run the Cleaning of the auto labeling
            if run_CADICA_3_AutoClean == 1
               CADICA_3_command = [CADICA_3_path, ' -c ', cadica_dir, ' -p ', subj_ids{j}, ' -s ', sess_ids{k}, ' -t ', task_names{l}, ' -i ', 'auto', ' -l ', num2str(lowfreq_cutoff), ' -h ', num2str(highfreq_cutoff), ' -b ', num2str(blur), ' < /dev/null'];
               fprintf(['Running: ', CADICA_3_command, '\n'])
               [status, cmdout_CADICA_3] = system(CADICA_3_command, '-echo');
               fprintf('Done with CADICA_3\n\n')
            end
            
            if run_CADICA_4_QC == 1
               fprintf('Running: CADICA_4_QC\n')
               prefix = ['sub-', subj_ids{j}, '_ses-', sess_ids{k}, '_task-', task_names{l}, '_'];
               suffix ='_bold.nii.gz';

               % Compare to 8p first (6 movement, CSF, and WM)
               CADICA_4_QC(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, prefix, suffix, '8p', 'auto', ic_filter_type, filterchoice, smoothchoice, lowfreq_cutoff, highfreq_cutoff, tr);
               
               % also compare to 9p (with global signal)
               CADICA_4_QC(cadica_dir, subj_ids{j}, sess_ids{k}, task_names{l}, prefix, suffix, '9p', 'auto', ic_filter_type, filterchoice, smoothchoice, lowfreq_cutoff, highfreq_cutoff, tr);
               fprintf('Done with CADICA_4\n\n')
            end
            fprintf(['Done with task ', task_names{l}, ' for session ', sess_ids{k}, ' for subject ', subj_ids{j}, '\n'])
            fprintf('\n')
            close all
            
        end
    end
end
fprintf('\n')


