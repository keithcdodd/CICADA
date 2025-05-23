function [Edge_r2, Edge_p, FD_r2, FD_p, DVARS_r2, DVARS_p, ...
    Outbrain_r2, Outbrain_p, WMCSF_r2, WMCSF_p, CSF_r2, CSF_p, Suscept_r2, Suscept_p, ...
    NotGM_r2, NotGM_p, GM_r2, GM_p] = get_r2_data(Edge_GM_corr, FD_GM_corr, ...
    DVARS_GM_corr, Outbrain_GM_corr, WMCSF_GM_corr, CSF_GM_corr, Suscept_GM_corr, NotGM_GM_corr, GM_GM_corr)

% Get R2 data from correlations
k = 2; % k is 2 for simple scatter or r from correlating two things

% calculate n for each sampling
Edge_n = length(Edge_GM_corr);
FD_n = length(FD_GM_corr);
DVARS_n = length(DVARS_GM_corr);
Outbrain_n = length(Outbrain_GM_corr);
WMCSF_n = length(WMCSF_GM_corr);
CSF_n = length(CSF_GM_corr);
Suscept_n = length(Suscept_GM_corr);
NotGM_n = length(NotGM_GM_corr);
GM_n = length(GM_GM_corr);

% Convert r to R2
Edge_r2 = Edge_GM_corr.^2;
FD_r2 = FD_GM_corr.^2;
DVARS_r2 = DVARS_GM_corr.^2;
Outbrain_r2 = Outbrain_GM_corr.^2;
WMCSF_r2 = WMCSF_GM_corr.^2;
CSF_r2 = CSF_GM_corr.^2;
Suscept_r2 = Suscept_GM_corr.^2;
NotGM_r2 = NotGM_GM_corr.^2;
GM_r2 = GM_GM_corr.^2;





% Calculate p value of mean r2 from beta distribution beta( (k-1)/2 , (n-k) / 2)
% good for comparing to null of no significant r2 (i.e., one sample
% testing). Is the mean R2 significant?
Edge_p = betacdf(mean(Edge_r2), (k-1)/2, (Edge_n-k)/2, 'upper');
FD_p = betacdf(mean(FD_r2), (k-1)/2, (FD_n-k)/2, 'upper');
DVARS_p = betacdf(mean(DVARS_r2), (k-1)/2, (DVARS_n-k)/2, 'upper');
Outbrain_p = betacdf(mean(Outbrain_r2), (k-1)/2, (Outbrain_n-k)/2, 'upper');
WMCSF_p = betacdf(mean(WMCSF_r2), (k-1)/2, (WMCSF_n-k)/2, 'upper');
CSF_p = betacdf(mean(CSF_r2), (k-1)/2, (CSF_n-k)/2, 'upper');
NotGM_p = betacdf(mean(NotGM_r2), (k-1)/2, (NotGM_n-k)/2, 'upper');
Suscept_p = betacdf(mean(Suscept_r2), (k-1)/2, (Suscept_n-k)/2, 'upper');
GM_p = betacdf(mean(GM_r2), (k-1)/2, (GM_n-k)/2, 'upper');


end