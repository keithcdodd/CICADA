function qc_corrs_table = grab_corr_sampling(sub_id, ses_id, task_name, Edge_Edge_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_Outbrain_corr, WMCSF_WMCSF_corr, CSF_CSF_corr, NotGM_NotGM_corr, GM_GM_corr, Suscept_Suscept_corr, samps)
% a function to grab samps number of random sampling from the region corrs.
% This makes it easier to plot group qc later (not overwhelming
% computationally)
initial_samps = samps;

% Make Group QC sampling reproducible by subject/session/task.
% Preserve the caller's random-number-generator state.
original_rng = rng;
rng_cleanup = onCleanup(@() rng(original_rng));

% Generate a stable seed from subject/session/task identifiers.
seed_key = char(string(sub_id) + "|" + ...
                string(ses_id) + "|" + string(task_name));

qc_seed = 5381;

for k = 1:numel(seed_key)
    qc_seed = mod(33 * qc_seed + double(seed_key(k)), 2^32);
end

rng(qc_seed, 'twister');
fprintf('Group QC sampling seed: %.0f\n', qc_seed);

% first create randperm of everything
NotGM_NotGM_randperm = randperm(size(NotGM_NotGM_corr, 1)); 
GM_GM_randperm = randperm(size(GM_GM_corr, 1));
DVARS_GM_randperm = randperm(size(DVARS_GM_corr, 1));
FD_GM_randperm = randperm(size(FD_GM_corr, 1));
CSF_CSF_randperm = randperm(size(CSF_CSF_corr, 1)); 
WMCSF_WMCSF_randperm = randperm(size(WMCSF_WMCSF_corr, 1));
Outbrain_Outbrain_randperm = randperm(size(Outbrain_Outbrain_corr, 1));
Edge_Edge_randperm = randperm(size(Edge_Edge_corr, 1)); 
Suscept_Suscept_randperm = randperm(size(Suscept_Suscept_corr, 1));

% Because we need everything to be the same size, we first test for
% smallest size or samps (which becomes our new samps):
samps = min([initial_samps, length(NotGM_NotGM_randperm), length(GM_GM_randperm), ...
    length(DVARS_GM_randperm), length(FD_GM_randperm), length(CSF_CSF_randperm), ...
    length(WMCSF_WMCSF_randperm), length(Outbrain_Outbrain_randperm), ...
    length(Edge_Edge_randperm), length(Suscept_Suscept_randperm)]);
% Now, we will sample at the same value for everyone, AND none will exceed
% the randperm limit! Ideally we sample the same for all group subjects,
% but this is close enough, very rarely will there not be enough samples,
% and if that does happen, it will be close

% now actually grab the sampling
NotGM_NotGM_randperm = NotGM_NotGM_randperm(1:samps);
GM_GM_randperm = GM_GM_randperm(1:samps);
DVARS_GM_randperm = DVARS_GM_randperm(1:samps);
FD_GM_randperm = FD_GM_randperm(1:samps);
CSF_CSF_randperm = CSF_CSF_randperm(1:samps);
WMCSF_WMCSF_randperm = WMCSF_WMCSF_randperm(1:samps);
Outbrain_Outbrain_randperm = Outbrain_Outbrain_randperm(1:samps);
Edge_Edge_randperm = Edge_Edge_randperm(1:samps);
Suscept_Suscept_randperm = Suscept_Suscept_randperm(1:samps);

NotGM_NotGM_corr = NotGM_NotGM_corr(NotGM_NotGM_randperm);
GM_GM_corr = GM_GM_corr(GM_GM_randperm);
DVARS_GM_corr = DVARS_GM_corr(DVARS_GM_randperm);
FD_GM_corr = FD_GM_corr(FD_GM_randperm);
CSF_CSF_corr = CSF_CSF_corr(CSF_CSF_randperm);
WMCSF_WMCSF_corr = WMCSF_WMCSF_corr(WMCSF_WMCSF_randperm);
Outbrain_Outbrain_corr = Outbrain_Outbrain_corr(Outbrain_Outbrain_randperm);
Edge_Edge_corr = Edge_Edge_corr(Edge_Edge_randperm);
Suscept_Suscept_corr = Suscept_Suscept_corr(Suscept_Suscept_randperm);

qc_corrs_array = {};
% Add to array!
for idx = 1:samps
    qc_corrs_array(idx,:) = {sub_id, ses_id, task_name, Edge_Edge_corr(idx), FD_GM_corr(idx), ...
        DVARS_GM_corr(idx), Outbrain_Outbrain_corr(idx), WMCSF_WMCSF_corr(idx), ...
        CSF_CSF_corr(idx), NotGM_NotGM_corr(idx), GM_GM_corr(idx), Suscept_Suscept_corr(idx)};
end

qc_corrs_labels = ["subject", "session", "task", "Edge_Edge_Corr", "FD_GM_Corr", ...
    "DVARS_GM_Corr", "Outbrain_Outbrain_Corr", "WMCSF_WMCSF_Corr", ...
    "CSF_CSF_Corr", "NotGM_NotGM_Corr", "GM_GM_Corr", "Suscept_Suscept_Corr"];
qc_corrs_table = array2table(qc_corrs_array, 'VariableNames',qc_corrs_labels);

end