function summary = cicada_build_support_products(task_dir, cicada_type, support_options)
% CICADA_BUILD_SUPPORT_PRODUCTS Build new products for an existing CICADA run.
%
% This does not rerun MELODIC, IC classification, or denoising.
if nargin < 3
    support_options = struct();
end
options = cicada_support_options(support_options);
if ~options.Enabled
    summary = table();
    return
end
if nargin < 2 || isempty(cicada_type)
    if isfile(fullfile(task_dir, 'ic_manual_selection', 'IC_manual_checker.csv'))
        cicada_type = 'manual';
    else
        cicada_type = 'auto';
    end
end

funcfile = fullfile(task_dir, 'funcfile.nii.gz');
funcmask = fullfile(task_dir, 'funcmask.nii.gz');
constrained = fullfile(task_dir, 'funcmask_constrained.nii.gz');
if ~isfile(constrained)
    constrained = make_constrained_funcmask(task_dir, funcfile, funcmask, 1, []);
end

make_funcdata_reliability(task_dir, funcfile, constrained, options);
products = cicada_locate_support_products(task_dir, cicada_type);
if isempty(products.SignalICOverlapRaw) || isempty(products.NoiseICOverlapRaw)
    error('CICADA:MissingRawComponentEvidence', ...
        'Could not locate SignalICOverlap and NoiseICOverlap for %s.', task_dir);
end
cicada_make_component_evidence(task_dir, cicada_type, ...
    products.SignalICOverlapRaw, products.NoiseICOverlapRaw, options);
summary = cicada_validate_support_products(task_dir, cicada_type);
end
