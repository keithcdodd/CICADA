function cicada_test_bids_metadata()
% Small synthetic, temporary fixture. Never touches real BIDS/CICADA data.
root = tempname;
mkdir(root);
cleanup = onCleanup(@() rmdir(root, 's')); %#ok<NASGU>
q = table({'001';'001';'002'}, {'01';'03';'03'}, ...
    {'rest';'rest';'rest'}, [2;5;7], ...
    'VariableNames', {'subject','session','task','number_total_ics'});

% No BIDS metadata: exact legacy columns and values.
[x, a] = cicada_group_bids_metadata(q, root);
assert(isequaln(x,q) && isempty(a.UsedFiles));

mkdir(fullfile(root,'sub-001'));
mkdir(fullfile(root,'sub-002'));
writeTSV(fullfile(root,'participants.tsv'), ...
    {'participant_id','group';'sub-001','treated';'sub-002','comparison'});
writeTSV(fullfile(root,'sub-001','sub-001_sessions.tsv'), ...
    {'session_id','treatment_condition';'ses-01','baseline';'ses-03','active'});
writeTSV(fullfile(root,'sub-002','sub-002_sessions.tsv'), ...
    {'session_id','treatment_condition';'ses-03','sham'});
mkdir(fullfile(root,'sub-002','ses-03'));
writeTSV(fullfile(root,'sub-002','ses-03', ...
    'sub-002_ses-03_scans.tsv'), ...
    {'filename','stimulation_on_during_scan'; ...
     'func/sub-002_ses-03_task-rest_bold.nii.gz','no'});
[y, b] = cicada_group_bids_metadata(q, root);
assert(isequal(y.number_total_ics, q.number_total_ics));
assert(isequal(string(y.bids_session_treatment_condition), ...
    ["baseline";"active";"sham"]));
assert(isequal(string(y.bids_participant_group), ...
    ["treated";"treated";"comparison"]));
assert(y.bids_scan_stimulation_on_during_scan(3) == "no");
assert(b.SessionMatches==3 && b.ScanMatches==1);

% Reject duplicate session IDs instead of assigning the wrong condition.
writeTSV(fullfile(root,'sub-002','sub-002_sessions.tsv'), ...
    {'session_id','treatment_condition';'ses-03','sham';'ses-03','active'});
failed = false;
try
    cicada_group_bids_metadata(q, root);
catch ME
    failed = strcmp(ME.identifier,'CICADA:BIDSMetadata:DuplicateKey');
end
assert(failed, 'Duplicate BIDS IDs must produce an error.');
fprintf('PASS: BIDS metadata absence, matches, scan-level fields, and duplicate-ID safety.\n');
end

function writeTSV(path, rows)
fid = fopen(path,'w');
assert(fid>=0);
for r = 1:size(rows,1)
    fprintf(fid,'%s\n',strjoin(rows(r,:),sprintf('\t')));
end
fclose(fid);
end
