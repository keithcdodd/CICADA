function [sample_times, info] = cicada_resolve_task_timing(TR, nVol, opts)
% CICADA_RESOLVE_TASK_TIMING  Resolve task-model sampling times safely.
%
%   [sample_times, info] = cicada_resolve_task_timing(TR, nVol)
%   [sample_times, info] = cicada_resolve_task_timing(TR, nVol, opts)
%
% Optional opts fields (all optional):
%   TaskSamplingReference - explicit scalar offset in [0,TR), or nVol-vector
%                           of monotonically increasing sample times (sec).
%   PreprocBoldJSON       - path to the JSON SIDECAR OF THE PREPROCESSED BOLD
%                           image, not the original BIDS acquisition JSON.
%   LegacyFMRIPrepSTC    - true ONLY when independent records confirm that
%                           older fMRIPrep applied slice-timing correction.
%   AcquisitionSliceTiming - original BIDS per-slice acquisition times (sec).
%                           Used only with LegacyFMRIPrepSTC=true, to build
%                           a PROVISIONAL midpoint reference estimate.
%
% Priority: explicit sampling reference > processed-BOLD metadata >
%           qualified legacy STC estimate > documented, unverified 0 s.
%
% IMPORTANT: A timing decision does NOT establish task-GLM inferential
% validity. A legacy midpoint reconstruction is only a candidate until
% software settings, BIDS event clock, and processing history are audited.
%
% This function does not read or write any imaging file and never alters
% the participant's data. Works without JSON for legacy datasets.

if nargin < 3 || isempty(opts)
    opts = struct();
end
assert(isstruct(opts) && isscalar(opts), ...
    'CICADA:TimingOptions', 'opts must be a scalar struct.');
assert(isnumeric(TR) && isscalar(TR) && isfinite(TR) && TR > 0, ...
    'CICADA:InvalidTR', 'TR must be a positive finite scalar.');
assert(isnumeric(nVol) && isscalar(nVol) && isfinite(nVol) && ...
    nVol >= 1 && nVol == round(nVol), ...
    'CICADA:InvalidNVol', 'nVol must be a positive integer.');

info = struct('Method', '', 'Confidence', '', ...
    'MetadataFile', '', 'OffsetSeconds', NaN, ...
    'ReferenceVerifiedFromMetadata', false, ...
    'Caveat', '', 'TRSeconds', double(TR), 'NumberOfVolumes', double(nVol));

% Highest priority: explicit user/pipeline time reference.
if isfield(opts, 'TaskSamplingReference') && ...
        ~isempty(opts.TaskSamplingReference)
    [sample_times, offset] = validate_reference( ...
        opts.TaskSamplingReference, TR, nVol);
    info.Method = 'explicit_reference';
    info.Confidence = 'user_or_pipeline_supplied_not_independently_verified';
    info.OffsetSeconds = offset;
    info.Caveat = ['The source of the explicit reference must be checked ', ...
        'against processed BOLD timing and the BIDS event onset clock.'];
    return;
end

% Optional modern fMRIPrep / derivative metadata.
if isfield(opts, 'PreprocBoldJSON') && ~isempty(opts.PreprocBoldJSON)
    jf = char(string(opts.PreprocBoldJSON));
    assert(isfile(jf), 'CICADA:TimingJSONMissing', ...
        'Specified preprocessed-BOLD JSON not found: %s', jf);
    meta = jsondecode(fileread(jf));
    assert(isstruct(meta) && isscalar(meta), 'CICADA:TimingJSONInvalid', ...
        'Preprocessed-BOLD JSON must contain an object.');
    info.MetadataFile = jf;

    if isfield(meta, 'RepetitionTime')
        metaTR = meta.RepetitionTime;
        assert(isnumeric(metaTR) && isscalar(metaTR) && ...
            isfinite(metaTR) && abs(metaTR-TR) <= max(1e-5, 1e-5*TR), ...
            'CICADA:TimingTRMismatch', ...
            'NIfTI TR does not match the preprocessed-BOLD JSON TR.');
    end

    hasSTC = isfield(meta, 'SliceTimingCorrected');
    hasStart = isfield(meta, 'StartTime');
    if hasSTC
        stc = meta.SliceTimingCorrected;
        assert((islogical(stc) || isnumeric(stc)) && isscalar(stc) && ...
            (isequal(stc, true) || isequal(stc, false)), ...
            'CICADA:TimingJSONInvalid', ...
            'SliceTimingCorrected must be a Boolean.');
        if logical(stc) && hasStart
            offset = meta.StartTime;
            assert(isnumeric(offset) && isscalar(offset) && ...
                isfinite(offset) && offset >= 0 && offset < TR, ...
                'CICADA:TimingJSONInvalid', ...
                'StartTime must be a finite offset in [0,TR).');
            sample_times = (0:nVol-1)' * TR + double(offset);
            info.Method = 'preproc_json_stc_starttime';
            info.Confidence = 'documented_processed_time_reference';
            info.OffsetSeconds = double(offset);
            info.ReferenceVerifiedFromMetadata = true;
            info.Caveat = ['StartTime gives the time reference declared by ', ...
                'the derivative. Event-clock alignment and GLM covariance ', ...
                'still require validation.'];
            return;
        end
        if ~logical(stc) && hasStart
            error('CICADA:TimingJSONContradiction', ...
                ['SliceTimingCorrected=false with StartTime present: ', ...
                 'do not infer an STC time reference automatically.']);
        end
        if logical(stc) && ~hasStart
            warning('CICADA:TimingSTCNoStartTime', ...
                ['Processed BOLD reports STC=true without StartTime. ', ...
                 'A usable common reference cannot be verified.']);
        end
        if ~logical(stc) && isfield(opts, 'LegacyFMRIPrepSTC') && ...
                isequal(opts.LegacyFMRIPrepSTC, true)
            error('CICADA:TimingJSONContradiction', ...
                ['Derivative reports STC=false but legacy options claim ', ...
                 'STC=true. Resolve conflicting provenance.']);
        end
        if ~logical(stc)
            warning('CICADA:TimingNoSTC', ...
                ['Processed BOLD reports STC=false. The nominal volume ', ...
                 'onset is a convention, not a common slice time.']);
        end
    else
        if hasStart
            warning('CICADA:TimingStartWithoutSTC', ...
                ['StartTime is present but correction status is absent. ', ...
                 'Not adopting an ambiguous reference automatically.']);
        else
            warning('CICADA:TimingJSONIncomplete', ...
                'JSON lacks the fields needed for an STC time reference.');
        end
    end
end

% Optional legacy fMRIPrep case. SliceTiming alone DOES NOT prove STC.
if isfield(opts, 'LegacyFMRIPrepSTC') && ...
        ~isempty(opts.LegacyFMRIPrepSTC) && ...
        isequal(opts.LegacyFMRIPrepSTC, true)
    if isfield(opts, 'AcquisitionSliceTiming') && ...
            ~isempty(opts.AcquisitionSliceTiming)
        st = opts.AcquisitionSliceTiming;
        assert(isnumeric(st) && isvector(st) && numel(st) >= 2 && ...
            all(isfinite(st(:))) && all(st(:) >= 0) && all(st(:) < TR), ...
            'CICADA:LegacySliceTimingInvalid', ...
            'AcquisitionSliceTiming must list valid per-slice times.');
        offset = (min(st(:)) + max(st(:))) / 2;
        sample_times = (0:nVol-1)' * TR + offset;
        info.Method = 'legacy_fmriprep_stc_midpoint_candidate';
        info.Confidence = 'provisional_reconstruction';
        info.OffsetSeconds = offset;
        info.Caveat = ['Legacy STC was reportedly verified. Midpoint ', ...
            'reconstruction is provisional until exact historical ', ...
            'reference behavior and run provenance are checked.'];
        warning('CICADA:LegacyMidpointEstimated', '%s', info.Caveat);
        return;
    else
        warning('CICADA:LegacySTCNoAcquisitionTimes', ...
            ['Legacy STC was declared, but slice times were not provided; ', ...
             'falling back to unverified nominal volume onset.']);
    end
end

% Nonblocking fallback for historical datasets with incomplete provenance.
sample_times = (0:nVol-1)' * TR;
info.Method = 'assumed_nominal_volume_onset';
info.Confidence = 'unverified';
info.OffsetSeconds = 0;
info.Caveat = ['No verified common sampling reference available. ', ...
    'Acceptable as a labeled screening convention; do not treat the ', ...
    'result as an inference-ready task GLM.'];
warning('CICADA:TaskTimingAssumed', '%s', info.Caveat);
end

function [times, offset] = validate_reference(reference, TR, nVol)
assert(isnumeric(reference) && isvector(reference) && ...
    ~isempty(reference) && all(isfinite(reference(:))), ...
    'CICADA:TimingReferenceInvalid', ...
    'TaskSamplingReference must be a finite scalar or vector.');
if isscalar(reference)
    offset = double(reference);
    assert(offset >= 0 && offset < TR, ...
        'CICADA:TimingReferenceInvalid', ...
        'Scalar offset must be in [0,TR).');
    times = (0:nVol-1)' * TR + offset;
else
    times = double(reference(:));
    assert(numel(times) == nVol && all(diff(times) > 0), ...
        'CICADA:TimingReferenceInvalid', ...
        'Sample times must have nVol entries and be increasing.');
    offset = NaN;
end
end
