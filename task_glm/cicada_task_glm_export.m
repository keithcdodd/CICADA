function P = cicada_task_glm_export(R,opts)
% CICADA_TASK_GLM_EXPORT  Optional export of an ALREADY-FITTED companion result.
%
% P = cicada_task_glm_export(R, opts)
%
% R must be the struct from cicada_task_glm with runAnalysis=true.
% opts.outputDir overrides the destination from the original preflight.
% DEFAULT opts.writeOutputs=false: preview ONLY; no files are written.
% To explicitly export, set opts.writeOutputs=true and select a NEW
% independent output directory. Sparse results additionally require
% opts.allowSparseTrial=true. Never overwrites existing output.
%
% This public entry point does not refit data, run FILM, or modify CICADA.
% Provisional conditional COPE/VARCOPE, NOT calibrated inference.

if nargin<2 || isempty(opts), opts=struct(); end
assert(isstruct(opts) && isscalar(opts), ...
    'opts must be a scalar struct');
assert(isstruct(R) && isscalar(R) && ...
    isfield(R,'fit') && isstruct(R.fit) && ...
    isfield(R,'preflight') && isstruct(R.preflight), ...
    'Supply a fitted result from cicada_task_glm(runAnalysis=true)');

if ~isfield(opts,'outputDir') || isempty(opts.outputDir)
    assert(isfield(R.preflight,'outputDir'), ...
        'Specify opts.outputDir');
    opts.outputDir=R.preflight.outputDir;
end
if ~isfield(opts,'writeOutputs') || isempty(opts.writeOutputs)
    opts.writeOutputs=false;
end
if ~isfield(opts,'allowSparseTrial') || isempty(opts.allowSparseTrial)
    opts.allowSparseTrial=false;
end

% Preserve validated writer implementation and its input/path guards.
helperDir=fullfile(fileparts(mfilename('fullpath')), 'helper_functions');
assert(isfolder(helperDir), ...
    'Missing task_glm/helper_functions. Keep companion folder intact');
originalPath=path;
restoreHelperPath=onCleanup(@()path(originalPath)); %#ok<NASGU>
addpath(helperDir,'-begin');
P=cicada_task_glm_reml_export_maps(R.fit,R.preflight,opts);
end
