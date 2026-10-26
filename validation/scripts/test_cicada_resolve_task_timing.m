function test_cicada_resolve_task_timing()
% TEST_CICADA_RESOLVE_TASK_TIMING  Local, nonimaging timing resolver tests.
% Run: test_cicada_resolve_task_timing
TR = 2; nVol = 172;
expected0 = (0:nVol-1)'*TR;

% The no-metadata path must work, but must never claim verification.
[t, p] = cicada_resolve_task_timing(TR, nVol);
assert(isequal(t, expected0));
assert(strcmp(p.Confidence, 'unverified'));
assert(~p.ReferenceVerifiedFromMetadata);

% Explicit scalar and equivalent vector must agree.
opts = struct('TaskSamplingReference', 0.97);
[ts, ps] = cicada_resolve_task_timing(TR, nVol, opts);
opts.TaskSamplingReference = expected0 + 0.97;
[tv, pv] = cicada_resolve_task_timing(TR, nVol, opts);
assert(isequal(ts, tv));
assert(strcmp(ps.Method, 'explicit_reference'));
assert(isnan(pv.OffsetSeconds));
assert(~ps.ReferenceVerifiedFromMetadata);

% Modern fMRIPrep derivative metadata with authoritative StartTime.
root = tempname;
mkdir(root);
cleanup = onCleanup(@() cleanup_test_dir(root)); %#ok<NASGU>
jsonfile = fullfile(root, 'sub-01_task-test_desc-preproc_bold.json');
meta = struct('RepetitionTime', TR, ...
    'SliceTimingCorrected', true, 'StartTime', 0.97);
write_json(jsonfile, meta);
[tj, pj] = cicada_resolve_task_timing(TR, nVol, ...
    struct('PreprocBoldJSON', jsonfile));
assert(isequal(tj, ts));
assert(pj.ReferenceVerifiedFromMetadata);
assert(strcmp(pj.Method, 'preproc_json_stc_starttime'));

% A user override must win even if a metadata file is also provided.
[tov, pov] = cicada_resolve_task_timing(TR, nVol, ...
    struct('PreprocBoldJSON', jsonfile, 'TaskSamplingReference', 0.1));
assert(isequal(tov, expected0+0.1));
assert(strcmp(pov.Method, 'explicit_reference'));

% No STC: not a verified shared slice reference.
meta = struct('RepetitionTime',TR,'SliceTimingCorrected',false);
write_json(jsonfile, meta);
[tno, pno] = cicada_resolve_task_timing(TR, nVol, ...
    struct('PreprocBoldJSON',jsonfile));
assert(isequal(tno, expected0));
assert(~pno.ReferenceVerifiedFromMetadata);

% Reported STC but missing StartTime -> labeled fallback.
meta = struct('RepetitionTime',TR,'SliceTimingCorrected',true);
write_json(jsonfile,meta);
[tmissing, pmissing] = cicada_resolve_task_timing(TR,nVol, ...
    struct('PreprocBoldJSON',jsonfile));
assert(isequal(tmissing,expected0));
assert(~pmissing.ReferenceVerifiedFromMetadata);

% Contradictory derivative metadata must be rejected.
meta = struct('RepetitionTime',TR, ...
    'SliceTimingCorrected',false,'StartTime',0.97);
write_json(jsonfile,meta);
assert_throws(@() cicada_resolve_task_timing(TR,nVol, ...
    struct('PreprocBoldJSON',jsonfile)), ...
    'CICADA:TimingJSONContradiction');

% TR mismatch must fail rather than silently shifting a design.
meta = struct('RepetitionTime',1.5, ...
    'SliceTimingCorrected',true,'StartTime',0.75);
write_json(jsonfile,meta);
assert_throws(@() cicada_resolve_task_timing(TR,nVol, ...
    struct('PreprocBoldJSON',jsonfile)), ...
    'CICADA:TimingTRMismatch');

% Legacy STC can be reconstructed only when independently declared.
legacy = struct('LegacyFMRIPrepSTC',true, ...
    'AcquisitionSliceTiming',[0,0.25,0.5,0.75,1.25,1.5]);
[tlegacy, pleg] = cicada_resolve_task_timing(TR,nVol,legacy);
assert(isequal(tlegacy,expected0+0.75));
assert(strcmp(pleg.Confidence,'provisional_reconstruction'));
assert(~pleg.ReferenceVerifiedFromMetadata);

% Acquisition SliceTiming ALONE must not be interpreted as STC proof.
[tonly, ponly] = cicada_resolve_task_timing(TR,nVol, ...
    struct('AcquisitionSliceTiming',legacy.AcquisitionSliceTiming));
assert(isequal(tonly,expected0));
assert(strcmp(ponly.Confidence,'unverified'));

fprintf('PASS: timing resolver tests (10 scenarios).\n');
end

function write_json(filename, meta)
fid=fopen(filename,'w');
assert(fid>0, 'Could not open temporary JSON fixture.');
closer=onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid,'%s',jsonencode(meta));
end

function assert_throws(f, identifier)
try
    f();
catch ME
    assert(strcmp(ME.identifier,identifier), ...
        'Expected error %s; received %s.',identifier,ME.identifier);
    return;
end
error('CICADA:MissingExpectedError', ...
    'Expected %s but no error occurred.',identifier);
end

function cleanup_test_dir(root)
if isfolder(root)
    rmdir(root,'s');
end
end
