function test_cicada_hrf_fractional()
% TEST_CICADA_HRF_FRACTIONAL Independent quadrature regression checks.
% Requires only cicada_task_hrf_design.m (plus standard MATLAB functions).
% No participant files, FSL, images, or writes.
TR = 2; nVol = 172;
p = [6 16 1 1 6 0 32];
names = {'short','long','early','impulse','late','outside'};
events = table( ...
    [6.3; 19.7; -3.25; 42.15; 330.3; 360], ...
    [1.1; 7.25; 6.2; 0; 30; 3], ...
    names', 'VariableNames', {'onset','duration','trial_type'});

X = cicada_task_hrf_design(events,TR,nVol,names,p,0);
assert(isequal(size(X),[172 6]));
assert(all(isfinite(X(:))));
assert(all(X(:,6)==0), 'Event entirely beyond run should not contribute.');
assert(any(X(:,3)~=0), 'Negative-onset event overlaps start of scan.');

% Reference constructed by numerical integration of the HRF PDF, not
% by the incomplete-gamma CDF used by CICADA.
Z = gammainc(p(7)/p(3),p(1)/p(3)) ...
    - gammainc(p(7)/p(4),p(2)/p(4))/p(5);
sample_times = (0:nVol-1)'*TR;
Xref = zeros(nVol,numel(names));
for j = 1:numel(names)
    for k = 1:nVol
        lag = sample_times(k) - events.onset(j) - p(6);
        d = events.duration(j);
        if d == 0
            Xref(k,j) = pdf_reference(lag,p,Z);
        else
            lo = max(0,lag-d);
            hi = min(p(7),lag);
            if hi > lo
                Xref(k,j) = integral( ...
                    @(u) pdf_reference(u,p,Z),lo,hi, ...
                    'AbsTol',1e-11,'RelTol',1e-10);
            end
        end
    end
end
err = max(abs(X(:)-Xref(:)));
assert(err < 1e-8, ...
    'Analytical continuous-time design differs from numerical quadrature.');

% Non-TR onset/duration should not be silently rounded.
Xshift = cicada_task_hrf_design(events,TR,nVol,names,p,0.91);
assert(any(abs(X(:)-Xshift(:))>1e-6));

fprintf('PASS: fractional HRF tests, max quadrature error %.3g.\n',err);
end

function h = pdf_reference(t,p,Z)
h = zeros(size(t));
active = t > 0 & t < p(7);
u = t(active);
a = p(1)/p(3); b = p(2)/p(4);
h(active) = (exp((a-1).*log(u)-u/p(3) ...
    - gammaln(a)-a*log(p(3))) ...
    - exp((b-1).*log(u)-u/p(4) ...
    - gammaln(b)-b*log(p(4)))/p(5))/Z;
end
