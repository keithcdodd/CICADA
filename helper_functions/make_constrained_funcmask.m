function funcmask_constrained = make_constrained_funcmask(output_dir, funcfile, funcmask, anatmask, use_kmeans, percent)
% This will make a constrained smaller funcmask then what is used in
% CICADA. The purpose is to use this constrained funcmask to help with
% group funcmask generation. This is needed because CICADA funcmask needs a
% large funcmask to help find sources of noise. However, this is not a
% great funcmask for analysis. This function makes a better funcmask for
% analysis.
% General idea is to remove darkest voxels from CICADA the tmean of the 
% funcfile  and then multiply by funcmask & anatmask. 
% NOTE: Anatmask should already be resampled to funcmask. e.g., with fsl flirt.
% Suggested to just use anatmask_resam in region_masks folder from CICADA
%
% Default method to find darkest voxels is with kmeans with 7 clusters.
% However, one can turn kmeans off (use_kmeans = 0;) and then give a percent
% (e.g. percent=20) to determine what percent cut off to decide are the
% darkest voxels ("susceptibility like") e.g., percent=20 would remove voxels below 20%
% quantile of funcfile that is masked by anatmask.
% use_kmeans and percent are both optional inputs. It is assumed use_kmeans is 1
% (do kmeans) and percent is not used ([])
% This can be helpful to run during standard CICADA, AND during group QC
% (so the user can easily change how the group_funcmask is calculated after
% CICADA is run).

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

if ~isfile(anatmask)
    fprintf(['Cannot find anatmask at ', anatmask, '\n'])
    return
end


cd(output_dir)

% Multiply funcfile by anatmask & then get temporal mean
funcfile_anatmask_command = ['fslmaths ', funcfile, ' -mul ', funcmask, ' -mul ', anatmask ' -Tmean ' output_dir, '/tmean_funcfile_anatmasked.nii.gz'];
call_fsl(funcfile_anatmask_command);
tmean_funcfile_anatmask = [output_dir, '/tmean_funcfile_anatmasked.nii.gz'];

% Read in file
tmean_funcfile_data = niftiread(tmean_funcfile_anatmask);

% grab data that is not masked out
tmean_vals = tmean_funcfile_data(tmean_funcfile_data > 0); % puts it into a list, only grabbing above 0

% do kmeans if percent value does not make sense
if ~exist('use_kmeans', 'var') && ~exist('percent', 'var')
    use_kmeans = 1; percent = [];
end

% only don't do kmeans if explicitely set to 0 as a double
if isa(use_kmeans, 'double') || use_kmeans ~= 0
    use_kmeans = 1;
end

% don't do kmeans, but no percent given, go to 20% default
if (use_kmeans == 0) && ~exist('percent', 'var')
    percent = 20; % default to 20% cut off, but this is kind of high
end

% basically, do kmeans unless kmeans==0 and percent value makes sense
if ~exist('percent', 'var') || use_kmeans ~= 0 || ~isa(percent, 'double') || percent > 99 || percent < 1
    fprintf('Doing kmeans to determine low data regions\n')
    % Get 5 points of quantiles
    Q = quantile(tmean_vals, [0, 0.17, 0.33, 0.5, 0.67, 0.83, 1]);
    % Now use kmeans clustering to find the darkest voxels (suscept-like):
    [val_idx, val_C] = kmeans(tmean_vals, 7, 'Start', [Q(1); Q(2); Q(3); Q(4); Q(5); Q(6); Q(7)], 'MaxIter', 10000);
    save('contrained_funcmask_vals', 'tmean_vals', 'Q', 'val_idx', 'val_C')
    min_voxel_value_idx = find(val_C == min(val_C));
    min_voxel_list_idx = val_idx == min_voxel_value_idx;
    voxel_val_cutoff = max(tmean_vals(min_voxel_list_idx)); % find the cut off value
else
    Q = quantile(tmean_vals, percent);
    voxel_val_cutoff = Q;

end

funcmask_constrained_data = single(tmean_funcfile_data > voxel_val_cutoff);

% The constrained funcmask will just be the opposite of min_voxels, you can
% write it out
funcmask_data_info = niftiinfo(funcmask);
niftiwrite(funcmask_constrained_data, 'funcmask_constrained', funcmask_data_info, 'Compressed', true)

funcmask_constrained = [output_dir, '/funcmask_constrained.nii.gz'];
