function F = cicada_task_glm_fingerprints(S)
% CICADA_TASK_GLM_FINGERPRINTS  Hash exactly the scientific inputs
% consumed by the Option 1 original-BOLD joint model, plus companion code.
% Read only; the cleaned BOLD image is not used by this GLM.
assert(isstruct(S) && isscalar(S) && isfield(S,'files'), ...
    'Supply a valid Option 1 preflight S');
F=struct();
F.schema='cicada_task_glm_fingerprints_v1';
F.companionVersion='0.3.3-dev-unreleased';
F.matlabVersion=version;
F.inputPaths=struct();
F.inputSHA256=struct();
inputs={'original','mask','mix','noise','rawHRF','cleanedJSON'};
for k=1:numel(inputs)
    key=inputs{k};
    assert(isfield(S.files,key),'Missing input field: %s',key);
    file=char(S.files.(key));
    F.inputPaths.(key)=char(java.io.File(file).getCanonicalPath());
    F.inputSHA256.(key)=cicada_task_glm_sha256(file);
end
hdir=fileparts(mfilename('fullpath'));
root=fileparts(hdir);
names={'cicada_task_glm', ...
    'cicada_task_glm_export', ...
    'cicada_task_glm_option1_film_pilot', ...
    'cicada_task_glm_reml_prepare', ...
    'cicada_task_glm_reml_fit_batch', ...
    'cicada_task_glm_reml_wholerun_readonly', ...
    'cicada_task_glm_reml_export_maps', ...
    'cicada_task_glm_sha256', ...
    'cicada_task_glm_fingerprints'};
F.codePaths=struct();
F.codeSHA256=struct();
for k=1:numel(names)
    key=names{k};
    if k<=2
        file=fullfile(root,[key '.m']);
    else
        file=fullfile(hdir,[key '.m']);
    end
    F.codePaths.(key)=char(java.io.File(file).getCanonicalPath());
    F.codeSHA256.(key)=cicada_task_glm_sha256(file);
end
end
