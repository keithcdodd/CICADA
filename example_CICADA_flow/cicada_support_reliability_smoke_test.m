function summary = cicada_support_reliability_smoke_test(task_dir, cicada_type, overwrite)
% CICADA_SUPPORT_RELIABILITY_SMOKE_TEST Build and validate rc1 products.
%
% summary = cicada_support_reliability_smoke_test(task_dir)
% summary = cicada_support_reliability_smoke_test(task_dir, 'auto', true)
%
% task_dir must contain an existing CICADA run. cicada_type is 'auto' or
% 'manual'. If omitted, the builder prefers manual when a manual checker is
% present. overwrite defaults to false.

if nargin < 2
    cicada_type = '';
end
if nargin < 3 || isempty(overwrite)
    overwrite = false;
end
validateattributes(overwrite, {'logical','numeric'}, {'scalar'});

options = struct();
options.Overwrite = logical(overwrite);
options.FailOnError = true;
summary = cicada_build_support_products(task_dir, cicada_type, options);

disp(summary)
out_dir = fullfile(task_dir, 'support_products');
if ~isfolder(out_dir)
    error('CICADA:SmokeTestNoOutput', 'support_products directory was not created.');
end
summary_path = fullfile(out_dir, 'support_reliability_smoke_test_summary.tsv');
writetable(summary, summary_path, 'FileType', 'text', 'Delimiter', '\t');
fprintf('Support/reliability smoke test passed. Summary: %s\n', summary_path);
end
