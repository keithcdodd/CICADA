% Run CADICA Cleaning Simpler script from Matlab (faster than calling in
% Bash)
clearvars

CADICAfol = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/derivatives/CADICA';
subjects = {'sub-102', 'sub-103', 'sub-105', 'sub-106', 'sub-107', 'sub-108', 'sub-109', ...
    'sub-112', 'sub-114', 'sub-115', 'sub-116', 'sub-117', 'sub-121', 'sub-122', 'sub-125', ...
    'sub-126', 'sub-128', 'sub-129', 'sub-130', 'sub-132', 'sub-134', 'sub-135', 'sub-136', ...
    'sub-138', 'sub-139', 'sub-140', 'sub-142', 'sub-143', 'sub-144', 'sub-145', 'sub-146', ...
    'sub-147', 'sub-149', 'sub-153', 'sub-156', 'sub-157', 'sub-158', 'sub-160', 'sub-161', ...
    'sub-164', 'sub-165', 'sub-168', 'sub-169', 'sub-171', 'sub-172', 'sub-173', 'sub-174', ...
    'sub-176', 'sub-178', 'sub-179', 'sub-181', 'sub-184', 'sub-185', 'sub-186', 'sub-187'};
sessions = {'ses-01', 'ses-02'};

subjectsconn = {'Subject001', 'Subject002', 'Subject003', 'Subject004', 'Subject005', 'Subject006', ...
    'Subject007', 'Subject008', 'Subject009', 'Subject010', 'Subject011', 'Subject012', 'Subject013', ...
    'Subject014', 'Subject015', 'Subject016', 'Subject017', 'Subject018', 'Subject019', 'Subject020', ...
    'Subject021', 'Subject022', 'Subject023', 'Subject024', 'Subject025', 'Subject026', 'Subject027', ...
    'Subject028', 'Subject029', 'Subject030', 'Subject031', 'Subject032', 'Subject033', 'Subject034', ...
    'Subject035', 'Subject036', 'Subject037', 'Subject038', 'Subject039', 'Subject040', 'Subject041', ...
    'Subject042', 'Subject043', 'Subject044', 'Subject045', 'Subject046', 'Subject047', 'Subject048', ...
    'Subject049', 'Subject050', 'Subject051', 'Subject052', 'Subject053', 'Subject054', 'Subject055'};
sessionsconn = {'Session001', 'Session002'};

connprojectdata = '/Volumes/VectoTec_VectoTech_Media_Rapid/AWESOME/Preproc_ICA_rest/conn_9P_restonly/RestNewICAImport/data';

cd(CADICAfol)

for j = 1:length(subjects)
    for k = 1:length(sessions)
        % go into correct folder
        cd(CADICAfol)
        cd(subjects{j})
        cd(sessions{k})
        
        clearvars -except CADICAfol subjects sessions j k connprojectdata subjectsconn sessionsconn
        % run script
        cadica_cleaning_simpler

        % Now reimport the covariates for the noise
        cd(connprojectdata)

        covariatefilename = ['COV_', subjectsconn{j}, '_', sessionsconn{k}, '.mat'];
        CADICAnoisefilename = [CADICAfol, '/', subjects{j}, '/', sessions{k}, '/CADICA_Noise.mat'];

        covariatefile = load(covariatefilename);
        CADICANoisefile = load(CADICAnoisefilename);
        index = find(contains(covariatefile.names, 'CADICA_Nois'));
        covariatefile.names{index} = 'CADICA_Noise';
        covariatefile.data{index} = CADICANoisefile.noise_covariates;
        save(covariatefilename, '-struct', 'covariatefile')


    end
end

cd(CADICAfol)


