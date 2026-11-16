function info = cicada_review_ics(taskDir, varargin)
% CICADA_REVIEW_ICS  Orthographic spatial review with native MELODIC plots.
%
% First launch:
%   cicada_review_ics(taskDir)
%   cicada_review_ics(taskDir, 'Anatomy', '/registered/anatomy.nii.gz')
%
% Primary classification workflow: edit IC_manual_checker.csv.
% Layout='auto' uses CICADA/templates/cicada_review.txt, when bundled;
% otherwise a personal FSLeyes layout or the default orthographic view.
% CICADA reads the bundled layout in FSLeyes itself, without installing it
% into the user's settings. 'Layout','melodic' still selects MELODIC mode.
% If exporting labels from FSLeyes, save to a new file before explicitly
% importing with 'ImportFSLeyesLabels' (this backs up the manual CSV).
%
% FSLeyes review shortcuts (when enabled): Ctrl+Shift+Right/Left or N/P to move
% forward/backward through the exact original CICADA PotentialICs order.
% These shortcuts are best-effort, dependent on FSLeyes/wxPython version.
%
% FSLeyes primarily displays maps and plots; IC_manual_checker.csv is the
% canonical manual decision file. Labels in FSLeyes are optional snapshots.
% ONLY the explicit ImportFSLeyesLabels action changes SignalLabel in
% IC_manual_checker.csv. No existing CICADA denoising code is changed.
% Never overwrites IC_auto_checker.csv. Creates a dated backup of an existing
% IC_manual_checker.csv before importing labels from FSLeyes.
%
% Default: bundled CICADA/templates/cicada_review.txt, when available;
% otherwise normal three-plane default or saved personal layout.
% The saved layout is applied as one FSLeyes operation, not built panel-by-panel.
% Uses session-level anatfol/anatfile.nii.gz when present; anatomical images
% must be aligned with the MELODIC IC image.
% FSLeyes auto-display handles positive and negative colour maps (~3-10
% initial range for IC maps with suitable amplitudes). This is visual only.
% MATLAB R2019b+ recommended. macOS/Linux; Windows launch unsupported.

if nargin < 1 || isempty(taskDir)
    taskDir = uigetdir(pwd, 'Select CICADA scan/task output folder');
    if isequal(taskDir, 0)
        info = struct('cancelled', true);
        return;
    end
end

taskDir = char(taskDir);
if ~isfolder(taskDir)
    error('CICADA:ReviewICs:TaskDir', 'Folder not found: %s', taskDir);
end

p = inputParser;
p.FunctionName = mfilename;
addParameter(p, 'Anatomy', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'OpenReport', false, @(x) islogical(x) && isscalar(x));
addParameter(p, 'OpenCSV', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'FSLeyesCommand', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ReviewShortcuts', true, @(x) islogical(x) && isscalar(x));
addParameter(p, 'Layout', 'auto', @(x) ischar(x) || (isstring(x) && isscalar(x)));
addParameter(p, 'ImportFSLeyesLabels', '', @(x) ischar(x) || (isstring(x) && isscalar(x)));
parse(p, varargin{:});
opts = p.Results;
layout = char(opts.Layout);
bundledLayout = '';
if strcmpi(layout, 'auto')
    % Source is CICADA/helper_functions/cicada_review_ics.m; the layout is
    % kept once in CICADA/templates, not copied into task folders.
    helperDir = fileparts(mfilename('fullpath'));
    templateFile = fullfile(fileparts(helperDir), 'templates', 'cicada_review.txt');
    if isfile(templateFile)
        % FSLeyes can apply a serialized layout directly from -r Python.
        % This also works when older FSLeyes versions use config.pkl.
        scene = 'default';
        bundledLayout = templateFile;
    else
        savedLayout = fullfile(getenv('HOME'), '.fsleyes', 'layouts', 'cicada_review.txt');
        if isfile(savedLayout)
            scene = 'cicada_review';
        else
            scene = 'default';  % No lightbox fallback.
        end
    end
elseif strcmpi(layout, 'ortho')
    scene = 'default';
else
    % Accept a user-chosen layout name, including built-in melodic.
    scene = layout;
end
if isempty(scene) || ~isempty(regexp(scene, '[^A-Za-z0-9_\-]', 'once'))
    error('CICADA:ReviewICs:Layout', ...
        'Layout name must contain only letters, digits, hyphens, or underscores.');
end

% Find the existing MELODIC folder without renaming or modifying it.
if isfolder(fullfile(taskDir, 'melodic'))
    melDir = fullfile(taskDir, 'melodic');
else
    melDir = taskDir;
    taskDir = fileparts(taskDir);
end

icPath = fullfile(melDir, 'melodic_IC.nii.gz');
if ~isfile(icPath)
    icPath = fullfile(melDir, 'melodic_IC.nii');
end
mixPath = fullfile(melDir, 'melodic_mix');
required = {icPath, mixPath, fullfile(melDir, 'melodic_FTmix')};
for i = 1:numel(required)
    if ~isfile(required{i})
        error('CICADA:ReviewICs:MissingMelodic', ...
            'Missing MELODIC file: %s', required{i});
    end
end

% IC_auto_checker.csv is already ordered by CICADA's ranking algorithm.
selectionDir = fullfile(taskDir, 'ic_auto_selection');
autoCsv = fullfile(selectionDir, 'IC_auto_checker.csv');
manualCsv = fullfile(selectionDir, 'IC_manual_checker.csv');
if ~isfile(autoCsv)
    error('CICADA:ReviewICs:MissingAuto', 'Not found: %s', autoCsv);
end

autoTable = readtable(autoCsv, 'TextType', 'string');
localCheckColumns(autoTable);
ids = localToNumbers(autoTable.PotentialICs);
autoLabels = localToNumbers(autoTable.SignalLabel);
nICs = size(readmatrix(mixPath), 2);
localValidateIDs(ids, autoLabels, nICs, autoCsv);

% Import is an EXPLICIT separate action to protect classifications.
importPath = char(opts.ImportFSLeyesLabels);
if ~isempty(importPath)
    [nChanged, backup] = localImport(importPath, manualCsv, autoTable, nICs);
    fprintf('Imported FSLeyes labels: %d changed vs. automatic CICADA.\n', nChanged);
    if ~isempty(backup)
        fprintf('Previous manual CSV backed up to: %s\n', backup);
    end
    fprintf('Updated manual CSV: %s\n', manualCsv);
    info = struct('cancelled', false, 'imported', true, ...
        'manualCsv', manualCsv, 'changedFromAuto', nChanged, 'backup', backup);
    return;
end

createdManual = false;
if ~isfile(manualCsv)
    [ok, msg] = copyfile(autoCsv, manualCsv);
    if ~ok
        error('CICADA:ReviewICs:Copy', 'Failed to create manual CSV: %s', msg);
    end
    createdManual = true;
end

% An existing manual CSV takes precedence so re-launching can resume work.
% Use autoTable for the ranking/evidence even if someone sorted manual CSV.
manualTable = readtable(manualCsv, 'TextType', 'string');
localCheckColumns(manualTable);
manualIDs = localToNumbers(manualTable.PotentialICs);
manualLabels = localToNumbers(manualTable.SignalLabel);
localValidateIDs(manualIDs, manualLabels, nICs, manualCsv);
[tf, location] = ismember(ids, manualIDs);
if ~all(tf)
    error('CICADA:ReviewICs:IDMismatch', ...
        'Manual CSV contains different component IDs from the automatic CSV.');
end
labelsInRankOrder = manualLabels(location);

reviewDir = fullfile(selectionDir, 'cicada_review');
if ~isfolder(reviewDir)
    mkdir(reviewDir);
end
queueCsv = fullfile(reviewDir, 'CICADA_review_order.csv');
labelsTxt = fullfile(reviewDir, 'CICADA_initial_fsleyes_labels.txt');
startupPy = fullfile(reviewDir, 'CICADA_fsleyes_startup.py');

% A human-readable reference queue, preserving all CICADA assessment columns.
queue = [table((1:nICs)', 'VariableNames', {'ReviewRank'}), autoTable];
queue.SignalLabel = labelsInRankOrder;
writetable(queue, queueCsv);

% FIX/Melview-compatible labels, indexed by ORIGINAL MELODIC IC number.
% These files are regenerated snapshots, never destinations for manual saves.
labelsByID = zeros(nICs, 1);
labelsByID(ids) = labelsInRankOrder;
localWriteFSLeyesLabels(labelsTxt, melDir, labelsByID);
[anatomy, anatomySource] = localFindAnatomy(taskDir, melDir, icPath, char(opts.Anatomy));
localWriteStartupScript(startupPy, labelsByID, ids, ...
    opts.ReviewShortcuts, anatomy, icPath, bundledLayout);

% FSLeyes launch: orthographic view by default; optional legacy MELODIC layout.
fsleyes = char(opts.FSLeyesCommand);
if isempty(fsleyes)
    [status, ~] = system('command -v fsleyes >/dev/null 2>&1');
    if status == 0
        fsleyes = 'fsleyes';
    else
        candidate = fullfile(getenv('FSLDIR'), 'bin', 'fsleyes');
        if isfile(candidate)
            fsleyes = candidate;
        end
    end
end

launched = false;
if isempty(fsleyes)
    warning('CICADA:ReviewICs:FSLeyes', ...
        'FSLeyes not found; pass ''FSLeyesCommand'' with its full path.');
elseif ispc
    warning('CICADA:ReviewICs:Windows', ...
        'FSLeyes launch on Windows is not implemented.');
else
    % FSLeyes loads this saved scene natively. Startup Python never changes
    % the view layout, so failed optional scripting cannot alter the panels.
    args = {localQuote(fsleyes), '-ad', '-s', scene, ...
        '-r', localQuote(startupPy)};
    if ~isempty(anatomy)
        args{end+1} = localQuote(anatomy);
    end
    args{end+1} = localQuote(melDir);
    logfile = fullfile(reviewDir, 'fsleyes_launch.log');
    cmd = [strjoin(args, ' '), ' >', localQuote(logfile), ' 2>&1 &'];
    launched = (system(cmd) == 0);
end

report = '';
candidates = {fullfile(melDir, 'report', '00index.html'), ...
    fullfile(melDir, 'report', 'index.html'), ...
    fullfile(melDir, 'report.html')};
for i = 1:numel(candidates)
    if isfile(candidates{i})
        report = candidates{i};
        break;
    end
end
if opts.OpenReport && ~isempty(report)
    localOpenFile(report);
end
if opts.OpenCSV && isfile(manualCsv)
    localOpenFile(manualCsv);
end

fprintf('\nCICADA IC review: %d components\n', nICs);
fprintf('  CICADA priority (first up to 15 ICs): %s\n', ...
    strjoin(arrayfun(@num2str, ids(1:min(15,end))', ...
    'UniformOutput', false), ', '));
fprintf('  Review sequence: %s\n', queueCsv);
fprintf('  Startup labels:  %s\n', labelsTxt);
fprintf('  Manual CSV:      %s\n', manualCsv);
fprintf('  FSLeyes scene:   %s\n', scene);
if ~isempty(bundledLayout)
    fprintf('  CICADA layout:   %s\n', bundledLayout);
elseif strcmp(scene, 'default')
    fprintf(['  To add the plots once, use View > Time series and\n' ...
        '  View > Power spectrum; save as layout cicada_review.\n']);
end
if isempty(anatomy)
    fprintf('  Background:      MELODIC mean EPI (no anatomical found)\n');
else
    fprintf('  Anatomy:         %s\n', anatomy);
    fprintf('  Anatomy source:  %s\n', anatomySource);
    fprintf('  IMPORTANT: Visually confirm registration to the IC map.\n');
end
if createdManual
    fprintf('  Created manual CSV as a copy of the automatic CSV.\n');
else
    fprintf('  Existing manual decisions were preserved and loaded.\n');
end
if opts.ReviewShortcuts
    fprintf('  FSLeyes Ctrl+Shift+Right/Left (or N/P): CICADA next/previous.\n');
    fprintf('  FSLeyes label snapshot is optional; edit the manual CSV as usual.\n');
else
    fprintf('  In FSLeyes, click Load labels and choose the startup labels file.\n');
end
fprintf(['  The manual CSV is the authoritative classification file.\n' ...
    '  FSLeyes Save labels/ImportFSLeyesLabels remains optional in melodic mode.\n' ...
    '  IC number is 1-based; FSLeyes volume index is 0-based.\n\n']);

info = struct('cancelled', false, 'imported', false, ...
    'taskDir', taskDir, 'melodicDir', melDir, 'autoCsv', autoCsv, ...
    'manualCsv', manualCsv, 'queueCsv', queueCsv, ...
    'labelsFile', labelsTxt, 'startupScript', startupPy, ...
    'reviewOrder', ids(:), 'report', report, ...
    'createdManualCsv', createdManual, 'fsleyesLaunchRequested', launched, ...
    'anatomy', anatomy, 'anatomySource', anatomySource, ...
    'layout', scene, 'bundledLayout', bundledLayout);
end

function localCheckColumns(t)
needed = {'PotentialICs', 'SignalLabel'};
if ~all(ismember(needed, t.Properties.VariableNames))
    error('CICADA:ReviewICs:Columns', ...
        'Expected columns PotentialICs and SignalLabel in the CICADA CSV.');
end
end

function nums = localToNumbers(values)
if isnumeric(values) || islogical(values)
    nums = double(values);
else
    nums = str2double(string(values));
end
nums = nums(:);
end

function localValidateIDs(ids, labels, nICs, path)
if numel(ids) ~= nICs || numel(labels) ~= nICs || ...
        any(~isfinite(ids)) || any(sort(ids(:)) ~= (1:nICs)') || ...
        any(~ismember(labels, [0, 1]))
    error('CICADA:ReviewICs:InvalidCSV', ...
        'Component IDs/SignalLabel must contain exactly 1:%d and binary labels: %s', ...
        nICs, path);
end
end

function localWriteFSLeyesLabels(filename, melDir, labelsByID)
fid = fopen(filename, 'w');
if fid < 0
    error('CICADA:ReviewICs:Write', 'Cannot create: %s', filename);
end
c = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', melDir);
bad = find(labelsByID == 0);
for i = 1:numel(labelsByID)
    if labelsByID(i) == 1
        fprintf(fid, '%d, Signal, False\n', i);
    else
        fprintf(fid, '%d, Unclassified noise, True\n', i);
    end
end
names = arrayfun(@num2str, bad, 'UniformOutput', false);
fprintf(fid, '[%s]\n', strjoin(names, ', '));
end

function localWriteStartupScript(filename, labelsByID, orderedIDs, ...
    reviewShortcuts, anatomy, icPath, bundledLayout)
% Apply the complete saved CICADA layout using FSLeyes' layout API;
% no manual panel-construction code or changes to FSLeyes user settings.
pyLabels = strjoin(arrayfun(@num2str, labelsByID(:)', ...
    'UniformOutput', false), ', ');
pyOrder = strjoin(arrayfun(@num2str, orderedIDs(:)', ...
    'UniformOutput', false), ', ');
lines = {
    'import wx'
    'import os'
    'import traceback'
    ['cicada_labels = [', pyLabels, ']']
    ['cicada_order = [', pyOrder, ']']
    ['cicada_shortcuts = ', localPyBool(reviewShortcuts)]
    ['cicada_anatomy = ', localPyString(anatomy)]
    ['cicada_ic_path = ', localPyString(icPath)]
    ['cicada_layout_path = ', localPyString(bundledLayout)]
    'cicada_state = {"index": 0, "tries": 0, "ready": False}'
    ''
    'if cicada_layout_path:'
    '    try:'
    '        from fsleyes import layouts as cicada_fsleyes_layouts'
    '        with open(cicada_layout_path, "r", encoding="utf-8") as f:'
    '            content = f.read()'
    '        title, separator, spec = content.partition("\n")'
    '        if not separator or "fsleyes.views.orthopanel.OrthoPanel" not in spec:'
    '            raise ValueError("Invalid bundled CICADA layout")'
    '        cicada_fsleyes_layouts.applyLayout(frame, "cicada_review", spec.strip())'
    '        print("CICADA: applied bundled layout: %s" % cicada_layout_path, flush=True)'
    '    except Exception as exc:'
    '        print("CICADA: bundled layout failed: %s" % exc, flush=True)'
    '        traceback.print_exc()'
    ''
    'def cicada_init():'
    '    if cicada_state["ready"]:'
    '        return'
    '    print("CICADA: open panels: %s" % [type(p).__name__ for p in frame.viewPanels], flush=True)'
    '    image = next((o for o in overlayList if getattr(o, "dataSource", None) and os.path.realpath(str(o.dataSource)) == os.path.realpath(cicada_ic_path)), None)'
    '    if image is None:'
    '        cicada_state["tries"] += 1'
    '        if cicada_state["tries"] < 15:'
    '            wx.CallLater(350, cicada_init)'
    '        else:'
    '            print("CICADA: IC overlay path not matched. Check log.", flush=True)'
    '            print("CICADA: loaded overlays: %s" % [str(getattr(o, "dataSource", None)) for o in overlayList], flush=True)'
    '        return'
    ''
    '    if cicada_anatomy:'
    '        try:'
    '            mean_path = image.getMeanFile()'
    '            for ov in overlayList:'
    '                if ov is not image and getattr(ov, "dataSource", None):'
    '                    if os.path.realpath(str(ov.dataSource)) == os.path.realpath(mean_path):'
    '                        displayCtx.getDisplay(ov).enabled = False'
    '        except Exception as exc:'
    '            print("CICADA: Could not hide MELODIC mean: %s" % exc, flush=True)'
    ''
    '    try:'
    '        volume_labels = overlayList.getData(image, "VolumeLabels")'
    '        if volume_labels is not None:'
    '            volume_labels.clear()'
    '            for idx, val in enumerate(cicada_labels):'
    '                volume_labels.addLabel(idx, "Signal" if val == 1 else "Unclassified noise")'
    '    except Exception as exc:'
    '        print("CICADA: optional FSLeyes labels unavailable: %s" % exc, flush=True)'
    ''
    '    displayCtx.selectOverlay(image)'
    '    def go_to_rank(rank):'
    '        if not 0 <= rank < len(cicada_order):'
    '            return'
    '        cicada_state["index"] = rank'
    '        ic = cicada_order[rank]'
    '        displayCtx.getOpts(image).volume = ic - 1'
    '        print("CICADA review rank %d/%d: IC %d" % (rank + 1, len(cicada_order), ic), flush=True)'
    ''
    '    if cicada_shortcuts:'
    '        def on_key(evt):'
    '            code = evt.GetKeyCode()'
    '            forward = code in (wx.WXK_RIGHT, ord("N"))'
    '            backward = code in (wx.WXK_LEFT, ord("P"))'
    '            if evt.ControlDown() and evt.ShiftDown() and (forward or backward):'
    '                go_to_rank(cicada_state["index"] + (1 if forward else -1))'
    '                return'
    '            evt.Skip()'
    '        frame.Bind(wx.EVT_CHAR_HOOK, on_key)'
    '        frame._cicada_review_keyhandler = on_key'
    ''
    '    cicada_state["ready"] = True'
    '    go_to_rank(0)'
    '    print("CICADA: Ready. Edit IC_manual_checker.csv; layout stays unchanged.", flush=True)'
    ''
    'cicada_init()'
};
fid = fopen(filename, 'w');
if fid < 0
    error('CICADA:ReviewICs:Write', 'Cannot create: %s', filename);
end
c = onCleanup(@() fclose(fid));
for i = 1:numel(lines)
    fprintf(fid, '%s\n', lines{i});
end
end

function out = localPyBool(value)
if value
    out = 'True';
else
    out = 'False';
end
end

function out = localPyString(value)
% Python raw string is unsafe for paths ending in a backslash; use repr-like
% single quote escaping, applied to all optional paths and filenames.
value = char(value);
value = strrep(value, '\', '\\');
value = strrep(value, char(39), [char(92), char(39)]);
out = [char(39), value, char(39)];
end

function [anatomy, source] = localFindAnatomy(taskDir, melDir, icPath, requested)
% Prioritize an explicitly supplied image. For the common MNI-preprocessed
% CICADA workflow, check the shared session-level anatfol first, then a
% matching-space TemplateFlow template; fall back to FSL standard if needed. Alignment MUST be visually
% checked before treating the background as an anatomical reference.
anatomy = '';
source = 'none';
if strcmpi(strtrim(requested), 'none')
    source = 'MELODIC mean EPI (requested)';
    return;
elseif ~isempty(requested)
    if ~isfile(requested)
        error('CICADA:ReviewICs:Anatomy', ...
            'Registered anatomical image not found: %s', requested);
    end
    anatomy = requested;
    source = 'user-specified';
    return;
end

% Prefer CICADA's shared anatomical image at the session level, without
% copying an anatomical NIfTI into each task output folder.
sessionAnat = fullfile(fileparts(taskDir), 'anatfol', 'anatfile.nii.gz');
if isfile(sessionAnat)
    anatomy = sessionAnat;
    source = 'CICADA session anatomical';
    return;
end

fslDir = getenv('FSLDIR');
tfDir = getenv('TEMPLATEFLOW_HOME');
homeDir = getenv('HOME');
locs = {fullfile(melDir, 'reg', 'standard.nii.gz')};
labels = {'MELODIC-registered standard'};
if ~isempty(tfDir)
    locs{end+1} = fullfile(tfDir, 'tpl-MNI152NLin2009cAsym', ...
        'tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz');
    labels{end+1} = 'TemplateFlow MNI152NLin2009cAsym';
end
if ~isempty(homeDir)
    locs{end+1} = fullfile(homeDir, '.cache', 'templateflow', ...
        'tpl-MNI152NLin2009cAsym', ...
        'tpl-MNI152NLin2009cAsym_res-02_T1w.nii.gz');
    labels{end+1} = 'TemplateFlow MNI152NLin2009cAsym';
end
if ~isempty(fslDir)
    locs{end+1} = fullfile(fslDir, 'data', 'standard', ...
        'MNI152_T1_2mm_brain.nii.gz');
    labels{end+1} = 'FSL MNI152 (approximate for MNI2009)';
end
for i = 1:numel(locs)
    if isfile(locs{i})
        anatomy = locs{i};
        source = labels{i};
        break;
    end
end
if isempty(anatomy)
    warning('CICADA:ReviewICs:NoAnatomy', ...
        ['No standard anatomical template was found. Using the MELODIC mean ' ...
         'EPI. Pass ''Anatomy'' with a registered T1 to view actual anatomy.']);
elseif ~strcmp(source, 'MELODIC-registered standard')
    warning('CICADA:ReviewICs:AnatomyAlignment', ...
        ['Selected standard-space template automatically; this is appropriate ' ...
         'ONLY when the IC image is registered to that template. Verify ' ...
         'alignment visually or pass a matching anatomical T1 via ''Anatomy''.']);
end
end

function [nChanged, backup] = localImport(filename, manualCsv, autoTable, nICs)
if ~isfile(filename)
    error('CICADA:ReviewICs:ImportNotFound', ...
        'FSLeyes export file not found: %s', filename);
end
[reviewLabels, seen] = localParseFSLeyesLabels(filename, nICs);
if ~all(seen)
    error('CICADA:ReviewICs:IncompleteImport', ...
        'The FSLeyes export must have an explicit label for every IC.');
end
ids = localToNumbers(autoTable.PotentialICs);
newLabels = reviewLabels(ids);
nChanged = sum(newLabels ~= localToNumbers(autoTable.SignalLabel));
updated = autoTable;
updated.SignalLabel = newLabels;

backup = '';
if isfile(manualCsv)
    backup = [manualCsv, '.backup_', datestr(now, 'yyyymmdd_HHMMSSFFF')];
    [ok, msg] = copyfile(manualCsv, backup);
    if ~ok
        error('CICADA:ReviewICs:Backup', ...
            'Could not back up existing manual CSV: %s', msg);
    end
end
% Write a temp file first, then atomically replace the working manual CSV.
tmp = [tempname(fileparts(manualCsv)), '.csv'];
writetable(updated, tmp);
[ok, msg] = movefile(tmp, manualCsv, 'f');
if ~ok
    error('CICADA:ReviewICs:ImportWrite', ...
        'Could not update manual CSV: %s', msg);
end
end

function [labels, seen] = localParseFSLeyesLabels(filename, nICs)
contents = fileread(filename);
rows = regexp(contents, '\r\n|\n|\r', 'split');
labels = nan(nICs, 1);
seen = false(nICs, 1);
badList = [];
for i = 1:numel(rows)
    line = strtrim(rows{i});
    if isempty(line)
        continue;
    end
    if startsWith(line, '[') && endsWith(line, ']')
        badText = strtrim(line(2:end-1));
        if ~isempty(badText)
            parts = strsplit(badText, ',');
            badList = str2double(strtrim(parts));
            if any(~isfinite(badList))
                error('CICADA:ReviewICs:BadList', ...
                    'FSLeyes bad-component list cannot be parsed.');
            end
        end
        continue;
    end
    match = regexp(line, '^(\d+)\s*,\s*(.+)$', 'tokens', 'once');
    if isempty(match)
        continue; % e.g., the MELODIC directory header
    end
    ic = str2double(match{1});
    if ~isfinite(ic) || ic < 1 || ic > nICs || seen(ic)
        error('CICADA:ReviewICs:ImportIDs', ...
            'Invalid/duplicate IC number in FSLeyes export: %s', line);
    end
    tokens = strtrim(strsplit(match{2}, ','));
    % Drop trailing classifier probability / noise boolean, if present.
    while ~isempty(tokens)
        last = strtrim(tokens{end});
        if any(strcmpi(last, {'True','False'})) || ...
                ~isnan(str2double(last))
            tokens(end) = [];
        else
            break;
        end
    end
    if isempty(tokens)
        error('CICADA:ReviewICs:ImportLabel', ...
            'No IC classification at component %d.', ic);
    end
    if any(strcmpi(tokens, 'Unknown'))
        error('CICADA:ReviewICs:Unknown', ...
            'Component %d is Unknown: classify it before importing.', ic);
    end
    labels(ic) = double(any(strcmpi(tokens, 'Signal')));
    seen(ic) = true;
end
if ~isempty(badList)
    expected = find(labels == 0)';
    if any(badList < 1) || any(badList > nICs) || ...
            numel(unique(badList)) ~= numel(badList) || ...
            ~isequal(sort(badList(:)), sort(expected(:)))
        error('CICADA:ReviewICs:BadListMismatch', ...
            'FSLeyes bad-component list conflicts with individual labels.');
    end
end
end

function quoted = localQuote(str)
q = char(39);
quoted = [q, strrep(char(str), q, [q, char(34), q, char(34), q]), q];
end

function localOpenFile(filename)
if ismac
    command = ['open ', localQuote(filename), ' >/dev/null 2>&1 &'];
elseif isunix
    command = ['xdg-open ', localQuote(filename), ' >/dev/null 2>&1 &'];
else
    try
        winopen(filename);
        return;
    catch
        edit(filename);
        return;
    end
end
if system(command) ~= 0
    warning('CICADA:ReviewICs:Open', 'Could not open: %s', filename);
end
end
