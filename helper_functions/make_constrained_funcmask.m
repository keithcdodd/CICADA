function funcmask_constrained = make_constrained_funcmask(output_dir, funcfile, funcmask, use_outliers, percent)
% This will make a constrained smaller funcmask then what is used in
% CICADA. The purpose is to use this constrained funcmask to help with
% group funcmask generation. This is needed because CICADA funcmask needs a
% large funcmask to help find sources of noise. However, this is not a
% great funcmask for analysis. This function makes a better funcmask for
% analysis.
% General idea is to remove darkest voxels from CICADA funcfile tmean (anatmasked)
%
% Default method to find darkest voxels is with matlab isoutlier.
% However, one can turn kmeans off (use_outliers = 0;) and then give a percent
% (e.g. percent=20) to determine what percent cut off to decide are the
% darkest voxels ("susceptibility like") e.g., percent=20 would remove voxels below 20%
% quantile of funcfile that is masked by anatmask.
% use_outliers and percent are both optional inputs. It is assumed use_outliers is 1
% (do isoutlier) and percent is not used ([])
% This can be helpful to run during standard CICADA, AND during group QC
% (so the user can easily change how the group_funcmask is calculated after
% CICADA is run).
% using isoutlier method will remove minimal amount of voxels

% get path to templates
helper_function_dir = fileparts(mfilename('fullpath')); % this gives current script path
cd([helper_function_dir, '/../'])
cicada_script_dir = pwd;
mni_anatmask = [cicada_script_dir, '/templates/mni_icbm152_nlin_asym_09c/mni_icbm152_t1_tal_nlin_asym_09c_mask.nii.gz'];

if ~isfolder(output_dir)
    mkdir(output_dir)
end

if ~isfile(funcfile)
    fprintf(['Cannot find funcfile at ', funcfile, '\n'])
    return
end

if ~isfile(funcmask)
    fprintf(['Cannot find funcmask at ', funcmask, '\n'])
    return
end

if ~isfile(mni_anatmask)
    fprintf(['Cannot find mni_anatmask at ', mni_anatmask, '\n'])
    return
end


cd(output_dir)

% first, flirt the mni_anatmask to funcfile space. mni_anatmask is better
% than subject anatmask in terms of retaining center slice of brain. Therefore, slightly better
% for group analysis.
command_1 = ['flirt -ref ', funcmask, ' -in ', mni_anatmask, ' -out ', output_dir, '/region_masks/mni_anatmask_resam.nii.gz -usesqform -applyxfm'];
command_2 = ['fslmaths ', output_dir, '/region_masks/mni_anatmask_resam.nii.gz -bin ', output_dir, '/region_masks/mni_anatmask_resam.nii.gz'];
call_fsl(command_1);
call_fsl(command_2);

% Multiply funcfile by mni_anatmask & then get temporal mean
funcfile_anatmask_command = ['fslmaths ', funcfile, ' -mul ', funcmask, ' -mul ', output_dir, '/region_masks/mni_anatmask_resam.nii.gz -Tmean ' output_dir, '/tmean_funcfile_anatmasked.nii.gz'];
call_fsl(funcfile_anatmask_command);
tmean_funcfile_anatmask = [output_dir, '/tmean_funcfile_anatmasked.nii.gz'];

% Read in file
tmean_funcfile_data = niftiread(tmean_funcfile_anatmask);

% grab data that is not masked out
tmean_vals = tmean_funcfile_data(tmean_funcfile_data > 0); % puts it into a list, only grabbing above 0

% do kmeans if percent value does not make sense
if ~exist('use_kmeans', 'var') && ~exist('percent', 'var')
    use_outliers = 1; percent = [];
end

% only don't do kmeans if explicitely set to 0 as a double
if isa(use_outliers, 'double') || use_outliers ~= 0
    use_outliers = 1;
end

% don't do kmeans, but no percent given, go to 20% default
if (use_outliers == 0) && ~exist('percent', 'var')
    percent = 20; % default to 20% cut off
end

% basically, do kmeans unless kmeans==0 and percent value makes sense
if ~exist('percent', 'var') || use_outliers ~= 0 || ~isa(percent, 'double') || percent > 99 || percent < 1
    fprintf('Determining low data regions. \n')
    % Use isoutlier to find minimum values within anatmasked funcfile
    low_signal_voxels = isoutlier(tmean_vals) & (tmean_vals < median(tmean_vals));
    if sum(low_signal_voxels == 0)
        voxel_val_cutoff = 0; % if there are no outliers, just keep it positive
    else
        voxel_val_cutoff = max(tmean_vals(low_signal_voxels));
    end
    
else
    Q = quantile(tmean_vals, percent/100);
    voxel_val_cutoff = Q;
end

funcmask_constrained_data = single(tmean_funcfile_data > voxel_val_cutoff);

% The constrained funcmask will just be the opposite of min_voxels, you can
% write it out
funcmask_data_info = niftiinfo(funcmask);
niftiwrite(funcmask_constrained_data, 'funcmask_constrained', funcmask_data_info, 'Compressed', true)

funcmask_constrained = [output_dir, '/funcmask_constrained.nii.gz'];

