function qc_corrs_table = grab_corr_sampling(Edge_Edge_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_Outbrain_corr, WMCSF_WMCSF_corr, CSF_CSF_corr, NotGM_NotGM_corr, GM_GM_corr, Suscept_Suscept_corr, samps)
% a function to grab samps number of random sampling from the region corrs.
% This makes it easier to plot group qc later (not overwhelming
% computationally)

NotGM_GM_randperm = randperm(size(NotGM_NotGM_corr, 1)); % and then grab ~samps of these
NotGM_GM_randperm = NotGM_GM_randperm(1:samps);

GM_GM_randperm = randperm(size(GM_GM_corr, 1)); % and then grab ~samps of these
GM_GM_randperm = GM_GM_randperm(1:samps);

DVARS_GMrandperm = randperm(size(DVARS_GM_corr, 1)); % and then grab ~samps of these
DVARS_GMrandperm = DVARS_GMrandperm(1:samps);

FD_GMrandperm = randperm(size(FD_GM_corr, 1)); % and then grab ~samps of these
FD_GMrandperm = FD_GMrandperm(1:samps);

CSF_GMrandperm = randperm(size(CSF_CSF_corr, 1)); % and then grab ~samps of these
CSF_GMrandperm = CSF_GMrandperm(1:samps);

WMCSF_GMrandperm = randperm(size(WMCSF_WMCSF_corr, 1)); % and then grab ~samps of these
WMCSF_GMrandperm = WMCSF_GMrandperm(1:samps);

Outbrain_GMrandperm = randperm(size(Outbrain_Outbrain_corr, 1)); % and then grab ~samps of these
Outbrain_GMrandperm = Outbrain_GMrandperm(1:samps);

Edge_GMrandperm = randperm(size(Edge_Edge_corr, 1)); % and then grab ~samps of these
Edge_GMrandperm = Edge_GMrandperm(1:samps);

Suscept_GMrandperm = randperm(size(Suscept_Suscept_corr, 1)); % and then grab ~samps of these
Suscept_GMrandperm = Suscept_GMrandperm(1:samps);

NotGM_NotGM_corr = NotGM_NotGM_corr(NotGM_GM_randperm);
GM_GM_corr = GM_GM_corr(GM_GM_randperm);
DVARS_GM_corr = DVARS_GM_corr(DVARS_GMrandperm);
FD_GM_corr = FD_GM_corr(FD_GMrandperm);
CSF_CSF_corr = CSF_CSF_corr(CSF_GMrandperm);
WMCSF_WMCSF_corr = WMCSF_WMCSF_corr(WMCSF_GMrandperm);
Outbrain_Outbrain_corr = Outbrain_Outbrain_corr(Outbrain_GMrandperm);
Edge_Edge_corr = Edge_Edge_corr(Edge_GMrandperm);
Suscept_Suscept_corr = Suscept_Suscept_corr(Edge_GMrandperm);

qc_corrs_array = [Edge_Edge_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_Outbrain_corr, WMCSF_WMCSF_corr, CSF_CSF_corr, NotGM_NotGM_corr, GM_GM_corr, Suscept_Suscept_corr];
qc_corrs_labels = {'Edge_Edge_Corr', 'FD_GM_Corr', ...
    'DVARS_GM_Corr', 'Outbrain_Outbrain_Corr', 'WMCSF_WMCSF_Corr', ...
    'CSF_CSF_Corr', 'NotGM_NotGM_Corr', 'GM_GM_Corr', 'Suscept_Suscept_Corr'};
qc_corrs_table = array2table(qc_corrs_array, 'VariableNames',qc_corrs_labels);

end