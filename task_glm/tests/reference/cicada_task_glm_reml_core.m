function R = cicada_task_glm_reml_core(Y, X0, C0, opts)
% CICADA_TASK_GLM_REML_CORE  Experimental original-time joint REML/FGLS.
%
% R = cicada_task_glm_reml_core(Y, X0, C0, opts)
%   Y:    N-by-V original CICADA input BOLD (NOT cleaned BOLD).
%   X0:   N-by-P full original-time joint design, INCLUDING an intercept,
%         unstandardized task HRFs, drift, and selected noise ICs.
%   C0:   K-by-P contrast matrix, with intercept weight in column 1.
%   opts.model: 'AR1' or 'AR2', explicitly required; NO model selection.
%
% Returns coefficients, COPE, *provisional* plug-in VARCOPE, REML grid
% estimates, scale and numerical diagnostics. Does not produce p-values,
% t/z maps, inferential claims, or any file I/O. Fixed-ICA assumptions only.
%
% This research core reproduces the covariance families and restricted
% profile score of the uploaded CICADA AR1/AR2 calibration utilities.
% Grid defaults: r1=-0.20:0.05:0.85; AR2 kappa2=-0.40:0.10:0.40.
% Test with up to 64 voxels by default. For future independent chunking,
% maxVoxels can be explicitly overridden after separate validation.
%
% EXAMPLE (after running the existing task GLM dry-run to obtain S.glm):
%   opts = struct('model','AR1');
%   R = cicada_task_glm_reml_core(Y64,S.glm.X0,S.glm.C0,opts);
%
% The AR1/AR2 covariance families are defined in the ORIGINAL volume-time
% coordinates. U below is a residual-contrast basis, NOT a time grid.

if nargin < 4 || isempty(opts), opts = struct(); end
assert(isstruct(opts) && isscalar(opts), 'opts must be a scalar struct');
assert(isfield(opts,'model') && ~isempty(opts.model), ...
    'Explicitly set opts.model to AR1 or AR2');
model = upper(char(string(opts.model)));
assert(ismember(model,{'AR1','AR2'}), 'Unknown covariance model');
opts = setDefault(opts,'r1Grid',-0.20:0.05:0.85);
opts = setDefault(opts,'k2Grid',-0.40:0.10:0.40);
opts = setDefault(opts,'maxVoxels',64);
assert(isscalar(opts.maxVoxels) && opts.maxVoxels >= 1 && ...
    (isinf(opts.maxVoxels) || ...
    (isfinite(opts.maxVoxels) && opts.maxVoxels==floor(opts.maxVoxels))), ...
    'Invalid maxVoxels');

Y = double(Y); X0 = double(X0); C0 = double(C0);
assert(ismatrix(Y) && ismatrix(X0) && ismatrix(C0) && ...
    all(isfinite(Y(:))) && all(isfinite(X0(:))) && ...
    all(isfinite(C0(:))), 'Inputs must be finite 2D numeric matrices');
[N,V] = size(Y);
[Nx,P] = size(X0);
assert(N >= 8 && V >= 1 && V <= opts.maxVoxels && Nx==N, ...
    'Y/X0 dimension mismatch or voxel count exceeds maxVoxels');
assert(P >= 2 && P < N && rank(X0)==P, ...
    'X0 must be full column rank with positive residual dimension');
assert(size(C0,2)==P && size(C0,1)>=1 && ...
    all(any(C0~=0,2)), 'Contrasts must be nonzero K-by-P rows');
% Our joint design convention includes an explicit constant column.
assert(max(abs(X0(:,1)-1))<1e-10, ...
    'Expected the intercept in the first column of X0');

rGrid = double(opts.r1Grid(:)');
if strcmp(model,'AR1')
    kGrid = 0;
else
    kGrid = double(opts.k2Grid(:)');
end
assert(~isempty(rGrid) && all(isfinite(rGrid)) && ...
    all(abs(rGrid)<1) && all(diff(rGrid)>0), ...
    'Invalid r1Grid');
assert(~isempty(kGrid) && all(isfinite(kGrid)) && ...
    all(abs(kGrid)<1) && all(diff(kGrid)>0), ...
    'Invalid k2Grid');

D = N-P;
U = null(X0');
assert(size(U,2)==D, 'Residual-contrast basis dimension mismatch');
Z = U'*Y;
assert(all(sum(Z.^2,1)>0), ...
    'Zero residual energy for at least one voxel');

% Design-aware Gaussian REML profile (up to theta-independent constants):
% log|U' A(theta) U| + D*log([z'(U'AU)^(-1)z]/D).
% Batch the V voxels: one covariance factorization per candidate grid pair.
bestScore = inf(1,V);
rHat = nan(1,V);
kHat = nan(1,V);
nCandidates = 0;
for r1 = rGrid
    for k2 = kGrid
        nCandidates = nCandidates+1;
        A = ar2Correlation(N,r1,k2);
        K = U'*A*U;
        K = (K+K')/2;
        [LK,flag] = chol(K,'lower');
        assert(flag==0, ...
            'Restricted covariance not positive definite at r1=%g k2=%g', ...
            r1,k2);
        W = LK\Z;
        q = sum(W.^2,1);
        assert(all(isfinite(q) & q>0), 'Invalid restricted quadratic form');
        score = 2*sum(log(diag(LK))) + D*log(q/D);
        take = score < bestScore;
        bestScore(take) = score(take);
        rHat(take) = r1;
        kHat(take) = k2;
    end
end
assert(all(isfinite(bestScore)) && all(isfinite(rHat)) && ...
    all(isfinite(kHat)), 'REML fitting failed');

% Cholesky GLS using ORIGINAL-time covariance. Whitening both X0 and Y
% means the estimator is comparable to the original-BOLD joint model.
Kc = size(C0,1);
Bhat = nan(P,V);
cope = nan(Kc,V);
varcope = nan(Kc,V);
sigma2 = nan(1,V);
scaledCond = nan(1,V);
rcondR = nan(1,V);
params = [rHat(:) kHat(:)];
[uniqueParams,~,groups] = unique(params,'rows');
for g = 1:size(uniqueParams,1)
    columns = find(groups==g)';
    A = ar2Correlation(N,uniqueParams(g,1),uniqueParams(g,2));
    LA = chol(A,'lower');
    Xw = LA\X0;
    Yw = LA\Y(:,columns);
    [~,RX] = qr(Xw,0);
    assert(rank(RX)==P, 'Whitened joint design is rank deficient');
    B = Xw\Yw;
    Ew = Yw-Xw*B;
    s2 = sum(Ew.^2,1)/D;
    assert(all(isfinite(s2) & s2>0), 'Invalid estimated residual scale');
    contrastGeometry = sum((RX'\C0').^2,1)';
    Bhat(:,columns) = B;
    cope(:,columns) = C0*B;
    varcope(:,columns) = contrastGeometry*s2;
    sigma2(columns) = s2;
    colScale = sqrt(sum(Xw.^2,1));
    xScaled = Xw./colScale;
    sv = svd(xScaled,'econ');
    scaledCond(columns) = sv(1)/sv(end);
    rcondR(columns) = rcond(RX);
end
assert(all(isfinite(Bhat(:))) && all(isfinite(cope(:))) && ...
    all(isfinite(varcope(:)) & varcope(:)>0), 'Invalid FGLS results');

R = struct();
R.status = 'PROVISIONAL_FIXED_ICA_CONDITIONAL';
R.model = model;
R.N = N;
R.nVoxels = V;
R.designRank = P;
R.residualDimension = D;
R.nContrasts = Kc;
R.nGridCandidates = nCandidates;
R.r1Grid = rGrid;
R.k2Grid = kGrid;
R.r1Hat = rHat;
R.kappa2Hat = kHat;
R.remlProfileScore = bestScore;
R.r1AtGridBoundary = (rHat==rGrid(1) | rHat==rGrid(end));
R.kappa2AtGridBoundary = false(1,V);
if strcmp(model,'AR2')
    R.kappa2AtGridBoundary = (kHat==kGrid(1) | kHat==kGrid(end));
end
R.beta = Bhat;
R.cope = cope;
R.varcopePlugIn = varcope;
R.residualScale2 = sigma2;
R.scaledWhitenedDesignCondition = scaledCond;
R.rcondWhitenedR = rcondR;
R.notes = ['Estimated covariance and fixed ICA labels are treated as ' ...
    'fixed for plug-in variance calculations. No t, p or z inference ' ...
    'is calibrated or authorized by this function.'];
end

function A = ar2Correlation(N,r1,k2)
% Stationary AR2 autocorrelation from lag-one ACF and lag-two PACF.
% k2=0 exactly gives AR1, as in the existing calibration scripts.
phi1 = r1*(1-k2);
phi2 = k2;
acf = zeros(N,1);
acf(1) = 1;
acf(2) = r1;
for t = 3:N
    acf(t) = phi1*acf(t-1)+phi2*acf(t-2);
end
A = toeplitz(acf);
end

function opts = setDefault(opts,name,value)
if ~isfield(opts,name) || isempty(opts.(name))
    opts.(name) = value;
end
end
