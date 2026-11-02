function cache = cicada_task_glm_reml_prepare(X0,C0,opts)
% CICADA_TASK_GLM_REML_PREPARE  Precompute one scan's fixed AR1/AR2 FGLS grid.
%
% cache = cicada_task_glm_reml_prepare(X0,C0,struct('model','AR1'));
% cache = cicada_task_glm_reml_prepare(X0,C0,struct('model','AR2'));
%
% X0 = N x P original-time joint design, intercept in column 1.
% C0 = K x P contrast matrix, with zero intercept contrast weight if desired.
% The design must include task HRFs, drift, and selected ICA nuisance EVs.
% The BOLD passed to the subsequent fit must be ORIGINAL CICADA input BOLD.
%
% This prepares covariance grid factorizations and GLS coefficient operators
% ONCE, for reuse across batches of voxels. No I/O, FILM, ICA, smoothing,
% automatic model selection, or statistical inference is performed.
% Plug-in VARCOPE from the fit remains PROVISIONAL (fixed-ICA assumptions).
%
% Defaults match cicada_task_glm_reml_core and earlier calibration scripts:
% r1 = -0.20:0.05:0.85, kappa2 = -0.40:0.10:0.40 (AR2 only).
%
% opts.maxBatchVoxels defaults to 64. It limits EACH subsequent fit call.
% opts.maxCacheBytes defaults to 256 MiB, a conservative cache guard.
%
% This is an additional experimental function; existing source is unchanged.

if nargin < 3 || isempty(opts), opts=struct(); end
assert(isstruct(opts) && isscalar(opts),'opts must be a scalar struct');
assert(isfield(opts,'model') && ~isempty(opts.model), ...
    'Specify opts.model as AR1 or AR2');
model=upper(char(string(opts.model)));
assert(ismember(model,{'AR1','AR2'}),'Unknown covariance model');
opts=setDefault(opts,'r1Grid',-0.20:0.05:0.85);
opts=setDefault(opts,'k2Grid',-0.40:0.10:0.40);
opts=setDefault(opts,'maxBatchVoxels',64);
opts=setDefault(opts,'maxCacheBytes',256*1024^2);
assert(isscalar(opts.maxBatchVoxels) && ...
    isfinite(opts.maxBatchVoxels) && opts.maxBatchVoxels>=1 && ...
    opts.maxBatchVoxels==floor(opts.maxBatchVoxels), ...
    'Invalid maxBatchVoxels');
assert(isscalar(opts.maxCacheBytes) && isfinite(opts.maxCacheBytes) ...
    && opts.maxCacheBytes>0,'Invalid maxCacheBytes');

X0=double(X0); C0=double(C0);
assert(isnumeric(X0) && isnumeric(C0) && ...
    ismatrix(X0) && ismatrix(C0) && ...
    all(isfinite(X0(:))) && all(isfinite(C0(:))), ...
    'X0 and C0 must be finite numeric matrices');
[N,P]=size(X0);
[nc,Pc]=size(C0);
assert(N>=8 && P>=2 && P<N && rank(X0)==P, ...
    'X0 must be N x P and full column rank with positive residual rank');
assert(Pc==P && nc>=1 && all(any(C0~=0,2)), ...
    'Contrasts must be nonzero K x P rows');
assert(max(abs(X0(:,1)-1))<1e-10, ...
    'Expected explicit intercept in X0(:,1)');

rGrid=double(opts.r1Grid(:)');
if strcmp(model,'AR1')
    kGrid=0;
else
    kGrid=double(opts.k2Grid(:)');
end
assert(~isempty(rGrid) && all(isfinite(rGrid)) && ...
    all(abs(rGrid)<1) && all(diff(rGrid)>0), ...
    'Invalid r1Grid');
assert(~isempty(kGrid) && all(isfinite(kGrid)) && ...
    all(abs(kGrid)<1) && all(diff(kGrid)>0), ...
    'Invalid k2Grid');
G=numel(rGrid)*numel(kGrid);
assert(G<=2000,'Too many covariance grid candidates');
D=N-P;
% Approximate persistent double-array size, excluding temporary factorizations.
cacheBytes=8*G*(D*D + P*N + nc + 4) + 8*(N*P + N*D + nc*P);
assert(cacheBytes<=opts.maxCacheBytes, ...
    'Cache requires about %.1f MiB, above configured maxCacheBytes', ...
    cacheBytes/1024^2);

U=null(X0');
assert(size(U,2)==D,'Residual subspace dimension mismatch');
cholRestricted=zeros(D,D,G);
betaOperator=zeros(P,N,G);
contrastGeometry=zeros(nc,G);
logDetRestricted=zeros(1,G);
scaledCondition=zeros(1,G);
rcondR=zeros(1,G);
params=zeros(G,2);

% Loop ordering matches the existing reference: outer r1, inner kappa2.
g=0;
for a=1:numel(rGrid)
    for b=1:numel(kGrid)
        g=g+1;
        r1=rGrid(a);
        k2=kGrid(b);
        params(g,:)=[r1 k2];
        A=ar2Correlation(N,r1,k2);
        [LA,flagA]=chol(A,'lower');
        assert(flagA==0, ...
            'Original-time covariance not positive definite at grid point %d',g);

        K=U'*A*U;
        K=(K+K')/2;
        [LK,flagK]=chol(K,'lower');
        assert(flagK==0, ...
            'Restricted covariance not positive definite at grid point %d',g);
        cholRestricted(:,:,g)=LK;
        logDetRestricted(g)=2*sum(log(diag(LK)));

        Xw=LA\X0;
        [Qw,Rw]=qr(Xw,0);
        assert(rank(Rw)==P, ...
            'Whitened design rank deficient at grid point %d',g);
        % Rw\Qw'\LA, evaluated without forming the inverse of A.
        % (LA'\Qw)' = Qw'/LA.
        betaOperator(:,:,g)=Rw\((LA'\Qw)');
        contrastGeometry(:,g)=sum((Rw'\C0').^2,1)';
        assert(all(contrastGeometry(:,g)>0), ...
            'Nonpositive contrast geometry at grid point %d',g);
        colNorm=sqrt(sum(Xw.^2,1));
        assert(all(colNorm>0),'Zero whitened design column');
        sv=svd(Xw./colNorm,'econ');
        scaledCondition(g)=sv(1)/sv(end);
        rcondR(g)=rcond(Rw);
    end
end

cache=struct();
cache.status='PREPARED_FIXED_ICA_CONDITIONAL';
cache.model=model;
cache.N=N;
cache.P=P;
cache.residualDimension=D;
cache.nContrasts=nc;
cache.nGridCandidates=G;
cache.r1Grid=rGrid;
cache.k2Grid=kGrid;
cache.params=params;
cache.C0=C0;
cache.U=U;
cache.cholRestricted=cholRestricted;
cache.logDetRestricted=logDetRestricted;
cache.betaOperator=betaOperator;
cache.contrastGeometry=contrastGeometry;
cache.scaledCondition=scaledCondition;
cache.rcondR=rcondR;
cache.maxBatchVoxels=opts.maxBatchVoxels;
cache.approxCacheBytes=cacheBytes;
cache.notes=['Read-only scan-level numerical cache. Design, original-BOLD ' ...
    'lineage, sample timing, ICA selection, and source provenance must be ' ...
    'verified outside this function. No inferential calibration is claimed.'];
end

function A=ar2Correlation(N,r1,k2)
% Same stationary AR2 correlation convention as the existing reference.
phi1=r1*(1-k2);
phi2=k2;
acf=zeros(N,1);
acf(1)=1;
acf(2)=r1;
for t=3:N
    acf(t)=phi1*acf(t-1)+phi2*acf(t-2);
end
A=toeplitz(acf);
end

function opts=setDefault(opts,name,value)
if ~isfield(opts,name) || isempty(opts.(name))
    opts.(name)=value;
end
end
