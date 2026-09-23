function [data_temporal, data_means, temporal_info] = ...
    cicada_temporal_transform(data, tr, fpass, detrend_degree, data_means)
% CICADA_TEMPORAL_TRANSFORM
%
% Apply CICADA temporal detrending / mean removal and optional Butterworth
% filtering to a matrix organized as:
%
%       [nTime x nSeries]
%
% IMPORTANT:
%   - Series means are REMOVED but are NOT added back here.
%   - This allows the identical temporal transform to be applied to BOLD
%     data and the complete MELODIC mixing matrix before matched ICA cleanup.
%   - A caller processing BOLD data may restore data_means afterward.
%
% 2F production policy (conservative hardening; validated algorithm kept):
%   - Missing/invalid detrend_degree defaults to degree 2, whether or not
%     temporal filtering is requested.
%   - Explicit detrend_degree <= 0 means mean removal only.
%   - If filtering is requested, polynomial detrending/mean removal occurs
%     FIRST, then the Butterworth filter is applied.
%   - Butterworth prototype order remains 2.
%   - Transfer-function [b,a] representation remains canonical.
%   - FILTFILT remains the zero-phase forward/backward implementation.
%   - No SOS substitution is made.
%
% Inputs:
%   data            nTime x nSeries numeric matrix
%   tr              repetition time in seconds
%   fpass           double [low_hz high_hz], or [] for no filtering
%                   [0 high] gives low-pass; [low Nyquist+] gives high-pass
%   detrend_degree  polynomial degree; default = 2 if missing/invalid
%   data_means      optional 1 x nSeries means supplied by caller
%
% Outputs:
%   data_temporal   temporally transformed data, mean not restored
%   data_means      1 x nSeries original means
%   temporal_info   provenance / validation structure

%% Inputs and basic validation

if nargin < 3
    fpass = [];
end

if nargin < 4
    detrend_degree = [];
end

if ~isnumeric(data) || ~ismatrix(data) || isempty(data)
    error('data must be a nonempty numeric [nTime x nSeries] matrix.');
end

if ~isreal(data) || any(~isfinite(data(:)))
    error('data must contain only finite real values.');
end

if ~isnumeric(tr) || ~isreal(tr) || ~isscalar(tr) || ...
        ~isfinite(tr) || tr <= 0
    error('tr must be a finite positive scalar in seconds.');
end

N = size(data, 1);
nSeries = size(data, 2);
fs = 1 / tr;
nyq = fs / 2;
run_duration_sec = N * tr;

if nargin < 5 || isempty(data_means)
    data_means = mean(data, 1);
else
    if ~isnumeric(data_means) || ~isreal(data_means) || ...
            numel(data_means) ~= nSeries || ...
            any(~isfinite(data_means(:)))
        error('data_means must contain one finite real mean per data series.');
    end
    data_means = reshape(data_means, 1, []);
end

% Preserve historical CICADA input contract: fpass is a MATLAB double pair.
% Malformed requests are treated as no filtering, but now this is explicit.
fpass_requested = ~isempty(fpass);
fpass_valid = isa(fpass, 'double') && isreal(fpass) && ...
    numel(fpass) == 2 && all(isfinite(fpass(:)));

if fpass_requested && ~fpass_valid
    warning('CICADA:InvalidFpass', ...
        ['fpass must be a finite real double [low_hz high_hz] pair. ' ...
         'Temporal filtering will be skipped.']);
    fpass = [];
end

filter_requested = fpass_valid && ~isempty(fpass);

%% Detrending / mean removal -- always before filtering

detrended = 0;
detrend_default_applied = false;

detrend_valid = isa(detrend_degree, 'double') && ...
    isreal(detrend_degree) && isscalar(detrend_degree) && ...
    isfinite(detrend_degree) && ~isempty(detrend_degree);

if ~detrend_valid
    detrend_default_applied = true;
    detrend_degree = 2;
    fprintf('   Detrending to default degree 2 polynomial...\n')
end

if detrend_degree > 0
    detrend_degree = round(detrend_degree);

    if ~detrend_default_applied
        fprintf('   Detrending to %d polynomial...\n', detrend_degree)
    end

    data_detrended = detrend(data, detrend_degree);
    detrended = 1;
    detrend_policy = 'polynomial';
else
    detrend_degree = 0;
    fprintf('   Mean removal only (no polynomial detrending).\n')
    data_detrended = data - data_means;
    detrended = 0;
    detrend_policy = 'mean_only';
end

%% Temporal filtering

bp = 0;
data_temporal = data_detrended;

requested_fpass = [];
effective_fpass = [];
filter_mode = 'none';
filter_skipped_reason = '';
cutoff_clamped_low = false;
cutoff_clamped_high = false;
filt_order = 2;
one_way_order = 0;
effective_forward_backward_order = 0;
filtfilt_minN_exclusive = NaN;
low_cutoff_period_sec = NaN;
low_cutoff_cycles_in_run = NaN;
short_run_low_cutoff_warning = false;

if filter_requested

    requested_fpass = reshape(double(fpass), 1, 2);
    low_hz  = requested_fpass(1);
    high_hz = requested_fpass(2);

    if low_hz < 0
        fprintf('   Low Hz cutoff is below 0; clamping to 0 Hz.\n')
        low_hz = 0;
        cutoff_clamped_low = true;
    end

    if high_hz > nyq
        fprintf(['   High Hz cutoff is above Nyquist; clamping to ', ...
            '%.12g Hz.\n'], nyq)
        high_hz = nyq;
        cutoff_clamped_high = true;
    end

    effective_fpass = [low_hz high_hz];

    % Harden impossible post-clamp ranges rather than silently creating an
    % unintended near-zero/Nyquist filter via coefficient clipping.
    if high_hz <= 0
        filter_skipped_reason = 'high_cutoff_not_above_zero';
        warning('CICADA:InvalidEffectiveFpass', ...
            'Effective high cutoff is <= 0 Hz. Filtering will be skipped.');

    elseif low_hz >= nyq
        filter_skipped_reason = 'low_cutoff_at_or_above_nyquist';
        warning('CICADA:InvalidEffectiveFpass', ...
            ['Effective low cutoff is at or above Nyquist. ' ...
             'Filtering will be skipped.']);

    else
        doHighPass = (low_hz > 0);
        doLowPass  = (high_hz < nyq);

        if doHighPass && doLowPass && (low_hz >= high_hz)
            filter_skipped_reason = 'low_cutoff_not_below_high_cutoff';
            warning('CICADA:InvalidEffectiveFpass', ...
                ['low_hz must be below high_hz for band-pass filtering. ' ...
                 'Filtering will be skipped.']);

        elseif ~doHighPass && ~doLowPass
            filter_skipped_reason = 'full_frequency_range';
            fprintf(['   Filter bounds imply the full sampled frequency ', ...
                'range; filtering will be skipped.\n']);

        else
            fprintf(['   Note: temporal filtering requires the MATLAB ', ...
                'Signal Processing Toolbox.\n'])

            if ~doHighPass && doLowPass
                filter_mode = 'lowpass';
                Wn = high_hz / nyq;
                Wn = min(max(Wn, 1e-6), 1 - 1e-6);
                [b, a] = butter(filt_order, Wn, 'low');
                one_way_order = filt_order;

            elseif doHighPass && ~doLowPass
                filter_mode = 'highpass';
                Wn = low_hz / nyq;
                Wn = min(max(Wn, 1e-6), 1 - 1e-6);
                [b, a] = butter(filt_order, Wn, 'high');
                one_way_order = filt_order;

            else
                filter_mode = 'bandpass';
                Wn = [low_hz high_hz] / nyq;
                Wn = min(max(Wn, 1e-6), 1 - 1e-6);
                [b, a] = butter(filt_order, Wn, 'bandpass');
                % MATLAB band-pass transformation doubles prototype order.
                one_way_order = 2 * filt_order;
            end

            % FILTFILT applies the transfer function forward and backward,
            % yielding a zero-phase response with squared magnitude.
            effective_forward_backward_order = 2 * one_way_order;

            filtfilt_minN_exclusive = ...
                3 * (max(length(a), length(b)) - 1) + 1;

            if N <= filtfilt_minN_exclusive
                error( ...
                    ['Time series too short for CICADA filtfilt ', ...
                     '(N=%d, need > %d samples for this filter).'], ...
                    N, filtfilt_minN_exclusive);
            end

            if doHighPass
                low_cutoff_period_sec = 1 / low_hz;
                low_cutoff_cycles_in_run = run_duration_sec * low_hz;

                if low_cutoff_cycles_in_run < 1
                    short_run_low_cutoff_warning = true;
                    warning('CICADA:ShortRunForLowCutoff', ...
                        ['Run duration is %.3f s, shorter than one period ', ...
                         '(%.3f s) at the %.6g Hz low cutoff. Filtering is ', ...
                         'still applied, but finite-run/edge behavior may ', ...
                         'strongly influence this cutoff.'], ...
                        run_duration_sec, low_cutoff_period_sec, low_hz);
                end
            end

            fprintf('  Filtering at %g %g Hz (%s)...\n', ...
                low_hz, high_hz, filter_mode)

            bp = 1;
            data_temporal = filtfilt(b, a, double(data_detrended));
        end
    end

else
    if fpass_requested && ~fpass_valid
        filter_skipped_reason = 'invalid_fpass';
    else
        filter_skipped_reason = 'not_requested';
    end
    fprintf('   Not applying temporal filtering.\n')
end

%% Return operation information / provenance

temporal_info = struct();

temporal_info.operation_order = 'detrend_then_filter';
temporal_info.detrended = detrended;
temporal_info.detrend_degree = detrend_degree;
temporal_info.detrend_policy = detrend_policy;
temporal_info.detrend_default_applied = detrend_default_applied;

temporal_info.filtered = bp;
temporal_info.filter_requested = filter_requested;
temporal_info.requested_fpass_hz = requested_fpass;
temporal_info.effective_fpass_hz = effective_fpass;
temporal_info.filter_mode = filter_mode;
temporal_info.filter_skipped_reason = filter_skipped_reason;
temporal_info.cutoff_clamped_low = cutoff_clamped_low;
temporal_info.cutoff_clamped_high = cutoff_clamped_high;

temporal_info.filter_family = 'Butterworth';
temporal_info.filter_representation = 'transfer_function_ba';
temporal_info.zero_phase_method = 'filtfilt';
temporal_info.prototype_order = filt_order;
temporal_info.one_way_filter_order = one_way_order;
temporal_info.effective_forward_backward_order = ...
    effective_forward_backward_order;

temporal_info.tr_sec = tr;
temporal_info.sampling_frequency_hz = fs;
temporal_info.nyquist_hz = nyq;
temporal_info.n_timepoints = N;
temporal_info.run_duration_sec = run_duration_sec;
temporal_info.filtfilt_minN_exclusive = filtfilt_minN_exclusive;
temporal_info.low_cutoff_period_sec = low_cutoff_period_sec;
temporal_info.low_cutoff_cycles_in_run = low_cutoff_cycles_in_run;
temporal_info.short_run_low_cutoff_warning = ...
    short_run_low_cutoff_warning;

end
