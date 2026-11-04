function [qc, audit] = cicada_group_bids_metadata(qc, bids_root)
% CICADA_GROUP_BIDS_METADATA  OPTIONAL BIDS metadata for Group CICADA QC.
%
% [qc, audit] = cicada_group_bids_metadata(qc, bids_root)
%
% READ-ONLY. Does not change existing QC values, subject/session selection,
% ICA inputs, or outlier decisions. If no matching standard metadata tables
% exist, the returned QC table is unchanged (same variables and values).
%
% Supported BIDS value tables (when present):
%   <root>/participants.tsv
%   <root>/sub-*/sub-*_sessions.tsv
%   <root>/sub-*/ses-*/sub-*_ses-*_scans.tsv
% For datasets without sessions, subject-level scans.tsv is also checked.
% The *_sessions.json and *_scans.json files (including root-level
% dictionaries) describe columns, not per-scan condition values; source
% dictionary paths are reported for provenance when available.
%
% Added columns use bids_participant_, bids_session_, bids_scan_ prefixes.
% Their VALUES are strings (including numeric-looking fields) so IDs,
% dates, and BIDS 'n/a' values are preserved during metadata attachment.
%
% Scan-level lookup is attempted only for unambiguously identified BOLD
% acquisitions by a Group CICADA task name such as 'rest' or
% 'foodpics_run-01'. Ambiguous/unmatched acquisitions remain unassigned.
%
% Validates unique BIDS participant_id, session_id, and scans filename.
% Does not implement arbitrary per-image JSON inheritance or use metadata
% to infer whether an acquisition is active/sham; treatment labels must
% be explicitly present in a TSV column.

if ~istable(qc)
    error('CICADA:BIDSMetadata:Input', 'First argument must be a QC table.');
end
required = {'subject','session','task'};
if ~all(ismember(required, qc.Properties.VariableNames))
    error('CICADA:BIDSMetadata:Keys', ...
        'QC table must contain subject, session, and task columns.');
end
if ~(ischar(bids_root) || (isstring(bids_root) && isscalar(bids_root)))
    error('CICADA:BIDSMetadata:Root', 'bids_root must be a path string.');
end
bids_root = char(bids_root);
if ~isfolder(bids_root)
    error('CICADA:BIDSMetadata:Root', 'BIDS root not found: %s', bids_root);
end
n = height(qc);
audit = struct('BIDSRoot', bids_root, 'UsedFiles', {{}}, ...
    'DictionaryFiles', {{}}, 'Warnings', {{}}, ...
    'ParticipantMatches', 0, 'SessionMatches', 0, ...
    'ScanMatches', 0, 'RowCount', n);
cache = containers.Map('KeyType', 'char', 'ValueType', 'any');

% Column-descriptor dictionaries contain column descriptions/levels, not
% the observation-specific values joined to Group CICADA QC rows.
dictNames = {'participants.json','sessions.json','scans.json'};
for d = 1:numel(dictNames)
    p = fullfile(bids_root, dictNames{d});
    if isfile(p)
        audit.DictionaryFiles{end+1} = p;
    end
end

participants = fullfile(bids_root, 'participants.tsv');
P = loadIfPresent(participants, 'participant_id');

for i = 1:n
    sub = "sub-" + erase(string(qc.subject(i)), "sub-");
    ses = "ses-" + erase(string(qc.session(i)), "ses-");
    task = string(qc.task(i));
    if ismissing(sub) || ismissing(ses) || ismissing(task) || ...
            strlength(sub) <= 4 || strlength(ses) <= 4
        audit.Warnings{end+1} = sprintf('Invalid QC identity at row %d', i);
        continue;
    end

    % Participant attributes: match exactly one row in participants.tsv.
    if ~isempty(P)
        ix = find(string(P.participant_id) == sub);
        if numel(ix) == 1
            attach(P, ix, 'participant_id', 'bids_participant_');
            audit.ParticipantMatches = audit.ParticipantMatches + 1;
            addUsed(participants);
        elseif isempty(ix)
            warnMissing(participants, char(sub));
        end
    end

    % Visit-level condition: distinct table per BIDS participant.
    sp = fullfile(bids_root, char(sub), sprintf('%s_sessions.tsv', char(sub)));
    S = loadIfPresent(sp, 'session_id');
    if ~isempty(S)
        ix = find(string(S.session_id) == ses);
        if numel(ix) == 1
            attach(S, ix, 'session_id', 'bids_session_');
            audit.SessionMatches = audit.SessionMatches + 1;
            addUsed(sp);
            localDict = fullfile(bids_root, char(sub), sprintf('%s_sessions.json', char(sub)));
            addDict(localDict);
        elseif isempty(ix)
            warnMissing(sp, char(ses));
        end
    end

    % Recording-level condition, if task/run maps to exactly one BIDS
    % _bold acquisition. Do not guess when several acquisitions match.
    scanFolder = fullfile(bids_root, char(sub), char(ses));
    scanFile = fullfile(scanFolder, sprintf('%s_%s_scans.tsv',char(sub),char(ses)));
    if ~isfile(scanFile)
        scanFolder = fullfile(bids_root, char(sub));
        scanFile = fullfile(scanFolder, sprintf('%s_scans.tsv',char(sub)));
    end
    C = loadIfPresent(scanFile, 'filename');
    if ~isempty(C)
        [ix, why] = matchBoldScan(string(C.filename), task);
        if numel(ix) == 1
            attach(C, ix, 'filename', 'bids_scan_');
            audit.ScanMatches = audit.ScanMatches + 1;
            addUsed(scanFile);
            addDict(fullfile(scanFolder, 'scans.json'));
            addDict(fullfile(scanFolder, sprintf('%s_%s_scans.json',char(sub),char(ses))));
            addDict(fullfile(scanFolder, sprintf('%s_scans.json',char(sub))));
        elseif ~isempty(why)
            audit.Warnings{end+1} = sprintf( ...
                '%s %s %s: %s (%s)', sub, ses, task, why, scanFile);
        end
    end
end

audit.UsedFiles = unique(audit.UsedFiles, 'stable');
audit.DictionaryFiles = unique(audit.DictionaryFiles, 'stable');
audit.Warnings = unique(audit.Warnings, 'stable');

    function T = loadIfPresent(p, key)
        p = char(p);
        if ~isfile(p)
            T = table();
            return;
        end
        if isKey(cache, p)
            T = cache(p);
            return;
        end
        % Read all metadata as text. Preserve leading zeros and 'n/a'.
        opts = detectImportOptions(p, 'FileType', 'text', ...
            'Delimiter', '\t', 'VariableNamingRule', 'preserve');
        opts = setvartype(opts, opts.VariableNames, 'string');
        T = readtable(p, opts);
        if ~ismember(key, T.Properties.VariableNames)
            error('CICADA:BIDSMetadata:KeyMissing', ...
                'Missing BIDS key %s in %s', key, p);
        end
        ids = string(T.(key));
        if any(ismissing(ids) | strlength(ids) == 0) || ...
                numel(unique(ids)) ~= numel(ids)
            error('CICADA:BIDSMetadata:DuplicateKey', ...
                'Missing or duplicated %s in %s', key, p);
        end
        cache(p) = T;
    end

    function attach(T, row, key, prefix)
        variables = T.Properties.VariableNames;
        for j = 1:numel(variables)
            name = variables{j};
            if strcmp(name, key)
                continue;
            end
            outName = matlab.lang.makeValidName([prefix name]);
            if ~ismember(outName, qc.Properties.VariableNames)
                qc.(outName) = repmat("n/a", n, 1);
            elseif ~isstring(qc.(outName))
                error('CICADA:BIDSMetadata:ColumnCollision', ...
                    'Existing QC column conflicts with BIDS metadata: %s',outName);
            end
            val = string(T{row, j});
            if isempty(val) || ismissing(val) || strlength(val) == 0
                val = "n/a";
            end
            qc.(outName)(i) = val;
        end
    end

    function addUsed(p)
        audit.UsedFiles{end+1} = char(p);
    end

    function addDict(p)
        if isfile(p)
            audit.DictionaryFiles{end+1} = char(p);
        end
    end

    function warnMissing(p, identifier)
        audit.Warnings{end+1} = sprintf('No %s entry in %s', identifier, p);
    end
end

function [idx, reason] = matchBoldScan(filename, task)
idx = [];
reason = '';
parts = regexp(char(task), '^([A-Za-z0-9]+)(?:_run-([A-Za-z0-9]+))?$', 'tokens', 'once');
if isempty(parts)
    reason = sprintf('Group task name cannot be safely resolved to a BIDS acquisition: %s',task);
    return;
end
label = string(parts{1});
runLabel = "";
if numel(parts) > 1
    runLabel = string(parts{2});
end
filename = replace(filename, '\', '/');
% Use actual _bold source filenames, not BIDS derivative/fMRIPrep files.
valid = endsWith(filename, '_bold.nii.gz') | endsWith(filename, '_bold.nii');
valid = valid & contains(filename, '_task-' + label + '_');
if strlength(runLabel) > 0
    valid = valid & contains(filename, '_run-' + runLabel + '_');
end
idx = find(valid);
if numel(idx) > 1
    idx = [];
    reason = 'multiple BIDS BOLD scans matched; need a unique task/run identifier';
end
end
