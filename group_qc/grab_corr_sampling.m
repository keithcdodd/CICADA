function qc_corrs_table = grab_corr_sampling(Edge_GM_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, NotGM_GM_corr, GM_GM_autocorr, samps)
% a function to grab samps number of random sampling from the region corrs.
% This makes it easier to plot group qc later (not overwhelming
% computationally)

NotGM_GM_randperm = randperm(size(NotGM_GM_corr, 1)); % and then grab ~samps of these
NotGM_GM_randperm = NotGM_GM_randperm(1:samps);

GM_GM_randperm = randperm(size(GM_GM_autocorr, 1)); % and then grab ~samps of these
GM_GM_randperm = GM_GM_randperm(1:samps);

DVARS_GMrandperm = randperm(size(DVARS_GM_corr, 1)); % and then grab ~samps of these
DVARS_GMrandperm = DVARS_GMrandperm(1:samps);

FD_GMrandperm = randperm(size(FD_GM_corr, 1)); % and then grab ~samps of these
FD_GMrandperm = FD_GMrandperm(1:samps);

CSF_GMrandperm = randperm(size(CSF_GM_corr, 1)); % and then grab ~samps of these
CSF_GMrandperm = CSF_GMrandperm(1:samps);

WMCSF_GMrandperm = randperm(size(WMCSF_GM_corr, 1)); % and then grab ~samps of these
WMCSF_GMrandperm = WMCSF_GMrandperm(1:samps);

Outbrain_GMrandperm = randperm(size(Outbrain_GM_corr, 1)); % and then grab ~samps of these
Outbrain_GMrandperm = Outbrain_GMrandperm(1:samps);

Edge_GMrandperm = randperm(size(Edge_GM_corr, 1)); % and then grab ~samps of these
Edge_GMrandperm = Edge_GMrandperm(1:samps);

NotGM_GM_corr = NotGM_GM_corr(NotGM_GM_randperm);
GM_GM_autocorr = GM_GM_autocorr(GM_GM_randperm);
DVARS_GM_corr = DVARS_GM_corr(DVARS_GMrandperm);
FD_GM_corr = FD_GM_corr(FD_GMrandperm);
CSF_GM_corr = CSF_GM_corr(CSF_GMrandperm);
WMCSF_GM_corr = WMCSF_GM_corr(WMCSF_GMrandperm);
Outbrain_GM_corr = Outbrain_GM_corr(Outbrain_GMrandperm);
Edge_GM_corr = Edge_GM_corr(Edge_GMrandperm);

qc_corrs_array = [Edge_GM_corr, FD_GM_corr, DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, NotGM_GM_corr, GM_GM_autocorr];
qc_corrs_labels = {'Edge_GM_Corr', 'FD_GM_Corr', ...
    'DVARS_GM_Corr', 'Outbrain_GM_Corr', 'WMCSF_GM_Corr', ...
    'CSF_GM_Corr', 'NotGM_GM_Corr', 'GM_GM_AutoCorr'};
qc_corrs_table = array2table(qc_corrs_array, 'VariableNames',qc_corrs_labels);

end