function X = cicada_task_hrf_design(events, TR, nVol, names, p, phase)
% Continuous-time double-gamma BIDS event design.
% Columns are unstandardized, amplitude-preserving regressors.
% Positive durations are boxcars; zero durations are impulses.


if nargin < 6 || isempty(phase)
    phase = 0;
    warning('CICADA:TimingAssumed', ...
        ['No sampling reference supplied. Using volume ', ...
         'onset times; timing is not independently verified.']);
end

assert(isfinite(TR) && TR > 0);
assert(nVol >= 1 && nVol == round(nVol));

if isscalar(phase)
    assert(isfinite(phase) && phase >= 0 && phase < TR);
    sample_times = (0:nVol-1)' * TR + phase;
else
    sample_times = double(phase(:));
    assert(numel(sample_times) == nVol);
    assert(all(isfinite(sample_times)));
    assert(all(diff(sample_times) > 0));
end


required = {'onset','duration','trial_type'};
assert(all(ismember(required, ...
    events.Properties.VariableNames)), ...
    'Events must contain onset, duration, trial_type.');

onset = double(events.onset);
duration = double(events.duration);
labels = strtrim(string(events.trial_type));
names = strtrim(string(names));

assert(all(isfinite(onset)) && all(isfinite(duration)));
assert(all(duration >= 0));

% Normalize the HRF to unit continuous-time area,
% truncated at p(7) seconds, as in the legacy model.
a1 = p(1)/p(3);
a2 = p(2)/p(4);
Z = gammainc(p(7)/p(3), a1) ...
    - gammainc(p(7)/p(4), a2)/p(5);
assert(Z > 0);

t = sample_times;
X = zeros(nVol, numel(names));

for k = 1:height(events)
    if strcmpi(strtrim(labels(k)), "baseline")
        continue
    end

    j = find(names == labels(k), 1);
    assert(~isempty(j), 'Unknown event condition.');

    lag = t - onset(k) - p(6);

    if duration(k) == 0
        response = hrf_pdf(lag, p, Z);
    else
        response = hrf_cdf(lag, p, Z) ...
            - hrf_cdf(lag-duration(k), p, Z);
    end

    X(:,j) = X(:,j) + response;
end
end

function F = hrf_cdf(t, p, Z)
u = max(0, min(t, p(7)));
F = (gammainc(u/p(3), p(1)/p(3)) ...
    - gammainc(u/p(4), p(2)/p(4))/p(5)) / Z;
end

function h = hrf_pdf(t, p, Z)
h = zeros(size(t));
idx = t > 0 & t < p(7);
u = t(idx);

a1 = p(1)/p(3);
a2 = p(2)/p(4);

g1 = exp((a1-1).*log(u) - u/p(3) ...
    - gammaln(a1) - a1*log(p(3)));
g2 = exp((a2-1).*log(u) - u/p(4) ...
    - gammaln(a2) - a2*log(p(4)));

h(idx) = (g1 - g2/p(5)) / Z;
end
