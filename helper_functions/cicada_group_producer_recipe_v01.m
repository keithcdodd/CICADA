function result = cicada_group_producer_recipe_v01(action, varargin)
% CICADA_GROUP_PRODUCER_RECIPE_V01
% External proposed production helper. No MRI modifications.
%
% Phase 1, immediately BEFORE matched Group preparation:
%   ticket = cicada_group_producer_recipe_v01('begin', cfg)
% Phase 2, immediately AFTER successful final Group BOLD preparation:
%   recipePath = cicada_group_producer_recipe_v01('finish', ...
%       ticket, preparedBOLD, actualTemporal, actualSpatial)
% Read-only verify later:
%   recipe = cicada_group_producer_recipe_v01('verify', recipePath)
%
% Requires cicada_task_glm_sha256 on MATLAB path. The producer-time
% handshake guarantees named input/source bytes did not change between
% begin and finish, and binds a resulting output file to those hashes.
% It does NOT prove numerical equivalence to the documented transform,
% historical lineage, validated inference, or protection from concurrent
% writers/power loss. Never use it to retroactively certify old Group files.
%
% cfg required fields: recipePath; preparationOutputDir; run (scalar struct); inputs (struct of
% named existing paths); codeFiles (struct of named existing paths);
% requestedProcessing (scalar struct). The input struct MUST contain:
% OriginalCopyBOLD, FullMelodicMix, NoiseSelection, FunctionalMask.

assert(exist('cicada_task_glm_sha256','file')==2, ...
    'The current CICADA task_glm/helper_functions path is required.');
assert(ischar(action)||(isstring(action)&&isscalar(action)), ...
    'Action must be scalar text.');
switch lower(char(action))
    case 'begin'
        assert(numel(varargin)==1, 'begin requires cfg.');
        cfg=varargin{1};
        assert(isstruct(cfg)&&isscalar(cfg), 'cfg must be a scalar struct.');
        required={'recipePath','preparationOutputDir','run','inputs','codeFiles','requestedProcessing'};
        for j=1:numel(required)
            assert(isfield(cfg,required{j}), 'Missing cfg.%s', required{j});
        end
        assert(isstruct(cfg.run)&&isscalar(cfg.run), 'Invalid run identity.');
        assert(isstruct(cfg.requestedProcessing)&&isscalar(cfg.requestedProcessing), ...
            'Invalid requested processing.');
        assert(isfield(cfg.requestedProcessing,'EffectiveProcessingPath') && ...
            strcmp(char(string(cfg.requestedProcessing.EffectiveProcessingPath)), ...
            'matched_temporal_ica'), ...
            'Only matched_temporal_ica is supported by this recipe helper.');
        p=char(string(cfg.recipePath));
        assert(~isempty(p) && ~isfile(p) && ~isfolder(p), ...
            'Recipe destination already exists or is invalid: %s',p);
        assert(isfolder(fileparts(p)), ...
            'Recipe parent directory must already exist: %s',fileparts(p));
        outDir=char(string(cfg.preparationOutputDir));
        assert(isfolder(outDir), 'Preparation output directory missing: %s',outDir);
        before=dir(outDir);
        initialNames={before(~ismember({before.name},{'.','..'})).name};
        expected={'OriginalCopyBOLD','FullMelodicMix','NoiseSelection','FunctionalMask'};
        assert(isstruct(cfg.inputs)&&isscalar(cfg.inputs), 'Invalid input map.');
        for j=1:numel(expected)
            assert(isfield(cfg.inputs,expected{j}), ...
                'Missing cfg.inputs.%s',expected{j});
        end
        result=struct();
        result.Schema='cicada_group_preparation_recipe_v01';
        result.Assurance='PRODUCER_TIME_HASH_LINK_NOT_NUMERICAL_RECONSTRUCTION';
        result.Run=cfg.run;
        result.RequestedProcessing=cfg.requestedProcessing;
        result.Inputs=local_capture(cfg.inputs);
        result.CodeFiles=local_capture(cfg.codeFiles);
        result.RecipePath=p;
        result.PreparationOutputDir=outDir;
        result.PreexistingOutputEntries=initialNames;
        result.StartedUTC=local_now_utc();
        result.MATLABVersion=version;

    case 'finish'
        assert(numel(varargin)==4, ...
            'finish requires ticket, preparedBOLD, actualTemporal, actualSpatial.');
        ticket=varargin{1};
        assert(isstruct(ticket)&&isscalar(ticket) && ...
            isfield(ticket,'Schema') && ...
            strcmp(ticket.Schema,'cicada_group_preparation_recipe_v01'), ...
            'Invalid producer ticket.');
        assert(isfield(ticket,'Assurance') && strcmp(ticket.Assurance, ...
            'PRODUCER_TIME_HASH_LINK_NOT_NUMERICAL_RECONSTRUCTION'), ...
            'Invalid ticket assurance.');
        preparedBOLD=char(string(varargin{2}));
        assert(isstruct(varargin{3}) && isscalar(varargin{3}), ...
            'Actual temporal metadata must be a scalar struct.');
        assert(isstruct(varargin{4}) && isscalar(varargin{4}), ...
            'Actual spatial metadata must be a scalar struct.');
        local_verify_entries(ticket.Inputs);
        local_verify_entries(ticket.CodeFiles);
        assert(strcmp(fileparts(preparedBOLD),char(ticket.PreparationOutputDir)), ...
            'Prepared output is outside the captured output directory.');
        [~,outputName,outputExt]=fileparts(preparedBOLD);
        assert(~ismember([outputName outputExt],ticket.PreexistingOutputEntries), ...
            'Prepared output existed before producer begin; refusing retroactive attestation.');
        assert(isfile(preparedBOLD), 'Prepared Group output is missing: %s', ...
            preparedBOLD);
        p=char(ticket.RecipePath);
        assert(~isfile(p)&&~isfolder(p), ...
            'Refusing to overwrite existing Group recipe: %s',p);
        result=ticket;
        result.PreparedBOLD=local_capture(struct('GroupPreparedCleanedBOLD',preparedBOLD));
        result.ActualTemporal=varargin{3};
        result.ActualSpatial=varargin{4};
        result.CompletedUTC=local_now_utc();
        result.Assurance='PRODUCER_TIME_HASH_LINK_NOT_NUMERICAL_RECONSTRUCTION';
        % New-file-only write; staging in the same directory as the recipe.
        stage=tempname(fileparts(p));
        cleanup=onCleanup(@() local_delete_stage(stage)); %#ok<NASGU>
        fid=fopen(stage,'w');
        assert(fid>=0,'Cannot open recipe staging file.');
        n=fprintf(fid,'%s\n',jsonencode(result));
        ok=fclose(fid);
        assert(n>0&&ok==0, 'Failed to write/close staged recipe.');
        assert(~isfile(p), 'Recipe destination appeared during staging.');
        [moved,msg]=movefile(stage,p);
        assert(moved,'Cannot finalize recipe: %s',msg);
        % Verify the actual serialized artifact and its linked files.
        cicada_group_producer_recipe_v01('verify',p);
        result=p;

    case 'verify'
        assert(numel(varargin)==1,'verify requires recipe path.');
        p=char(string(varargin{1}));
        assert(isfile(p),'Recipe file missing: %s',p);
        result=jsondecode(fileread(p));
        assert(isfield(result,'Schema') && strcmp(result.Schema, ...
            'cicada_group_preparation_recipe_v01'), ...
            'Refusing non-producer or unknown recipe schema.');
        assert(isfield(result,'Assurance') && strcmp(result.Assurance, ...
            'PRODUCER_TIME_HASH_LINK_NOT_NUMERICAL_RECONSTRUCTION'), ...
            'Refusing recipe with wrong assurance status.');
        assert(isfield(result,'StartedUTC') && ...
            isfield(result,'CompletedUTC') && ...
            isfield(result,'ActualTemporal') && ...
            isfield(result,'ActualSpatial'), ...
            'Incomplete producer recipe.');
        local_verify_entries(result.Inputs);
        local_verify_entries(result.CodeFiles);
        local_verify_entries(result.PreparedBOLD);
        fprintf('PASS: producer-time input/output hash link: %s\n',p);
        fprintf(['LIMIT: no numerical reconstruction or inferential ' ...
            'calibration is established.\n']);
    otherwise
        error('Unknown action: %s',char(action));
end
end

function entries=local_capture(paths)
assert(isstruct(paths)&&isscalar(paths) && ~isempty(fieldnames(paths)), ...
    'File map must be a nonempty scalar struct.');
keys=fieldnames(paths);
entries=struct();
for k=1:numel(keys)
    path=char(string(paths.(keys{k})));
    assert(isfile(path),'Required producer file is missing: %s',path);
    d=dir(path);
    entries.(keys{k})=struct('path',path,'sizeBytes',double(d.bytes), ...
        'sha256',cicada_task_glm_sha256(path));
end
end

function local_verify_entries(entries)
assert(isstruct(entries)&&isscalar(entries)&&~isempty(fieldnames(entries)), ...
    'Malformed recipe file entries.');
keys=fieldnames(entries);
for k=1:numel(keys)
    e=entries.(keys{k});
    assert(isfield(e,'path')&&isfield(e,'sizeBytes')&&isfield(e,'sha256'), ...
        'Incomplete recipe file entry: %s',keys{k});
    p=char(string(e.path));
    assert(isfile(p),'Required producer file now missing: %s',p);
    d=dir(p);
    assert(double(d.bytes)==double(e.sizeBytes), ...
        'Producer file size changed: %s',p);
    assert(strcmpi(char(e.sha256),cicada_task_glm_sha256(p)), ...
        'Producer file SHA-256 changed: %s',p);
end
end

function local_delete_stage(p)
if isfile(p), delete(p); end
end

function s=local_now_utc()
t=datetime('now','TimeZone','UTC');
t.Format="yyyy-MM-dd'T'HH:mm:ss";
s=[char(t) 'Z'];
end
