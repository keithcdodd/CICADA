% CADICA_4_QC
% CADICA_3_Clean calculated parcellation text files from Gordon at last of
% denoised files. This script just reads them and runs correlations to plot
% QC
home='/home/keithdodd/ExampleDataLocal/CADICA_Updated';
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
        cd 'QC/parctimeseries'
        
        % Get GM prob mask and create GM mask and not GM mask

        % Select the timeseries for GM and not GM

        % randperm selection for both

        % corr top 100 for both

        % corr just GM to itself

        % grab lower triangles and plot them

        % GM to GM should skew left, GM to out of GM should center 0


        test = load('ts_1parc_s_bp_ICADenoised.txt');
        ts_ICADenoised = zeros(length(test), 300); % using schaefer 300 parcellations
        ts_9p = ts_ICADenoised;
        ts_4r = ts_ICADenoised;

        for j = 1:300
            currfile_ICADenoised = ['ts_', num2str(j), 'parc_s_bp_ICADenoised.txt'];
            currfile_9p = ['ts_', num2str(j), 'parc_s_bp_9p.txt'];
            currfile_4r = ['ts_', num2str(j), 'parc_s_bp_4r.txt'];
            ts_ICADenoised(:,j) = load(currfile_ICADenoised);
            ts_9p(:,j) = load(currfile_9p);
            ts_4r(:,j) = load(currfile_4r);
        end

        corr_mat_ICADenoised = corr(ts_ICADenoised,ts_ICADenoised);
        low_tri_ICADenoised = tril(corr_mat_ICADenoised,-1);

        corr_mat_9p = corr(ts_9p,ts_9p);
        low_tri_9p = tril(corr_mat_9p,-1);

        corr_mat_4r = corr(ts_4r,ts_4r);
        low_tri_4r = tril(corr_mat_4r,-1);

    end
end

figure
histogram(low_tri_ICADenoised(low_tri_ICADenoised ~= 0))

figure
histogram(low_tri_9p(low_tri_9p ~= 0))

figure
histogram(low_tri_4r(low_tri_4r ~= 0))